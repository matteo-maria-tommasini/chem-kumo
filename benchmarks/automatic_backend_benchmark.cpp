#include "chem/kumo.hpp"
#include "automatic_backend.hpp"
#include "evaluation_plan.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;
namespace {
using Clock = std::chrono::steady_clock;
volatile double sink = 0.0;
constexpr std::array components{density::total, density::alpha, density::beta, density::spin};

enum class Workload { Density, ELF, Hessian };
const char* name(Workload w) {
    switch (w) {
        case Workload::Density: return "density";
        case Workload::ELF: return "elf_ingredients";
        case Workload::Hessian: return "hessian_tau";
    }
    return "unknown";
}
Request request_for(Workload w) {
    Request req;
    for (auto d : components) {
        req.add(d.value());
        if (w != Workload::Density) req.add(d.gradient()).add(d.tau());
        if (w == Workload::Hessian) req.add(d.hessian()).add(d.laplacian());
    }
    return req;
}
const char* backend_name(Backend b) {
    switch (b) {
        case Backend::AUTO: return "AUTO";
        case Backend::AO_P: return "AO_P";
        case Backend::AO_MO: return "AO_MO";
    }
    return "UNKNOWN";
}
struct Comparison { double max_abs = 0, max_scaled = 0, checksum_ao = 0, checksum_mo = 0; };
Comparison compare(const Results& ao, const Results& mo, const Request& req) {
    Comparison result;
    for (const auto key : req.properties()) {
        const auto& a = ao.get(key);
        const auto& b = mo.get(key);
        if (a.points != b.points || a.components.size() != b.components.size())
            throw std::runtime_error("different result shapes");
        for (std::size_t c = 0; c < a.components.size(); ++c) {
            const auto av = a.component(c), bv = b.component(c);
            for (std::size_t p = 0; p < a.points; ++p) {
                if (!std::isfinite(av[p]) || !std::isfinite(bv[p]))
                    throw std::runtime_error("non-finite field value");
                const double diff = std::abs(av[p] - bv[p]);
                result.max_abs = std::max(result.max_abs, diff);
                result.max_scaled = std::max(result.max_scaled,
                    diff / std::max({1.0, std::abs(av[p]), std::abs(bv[p])}));
                result.checksum_ao += av[p];
                result.checksum_mo += bv[p];
            }
        }
    }
    return result;
}
double checksum(const Results& result, const Request& req) {
    double sum = 0.0;
    for (const auto key : req.properties()) {
        const auto& field = result.get(key);
        for (const auto& plane : field.components)
            for (double value : plane) sum += value;
    }
    if (!std::isfinite(sum)) throw std::runtime_error("non-finite checksum");
    return sum;
}
template<class F> double median_seconds(F&& f, int repeats) {
    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(repeats));
    for (int i = 0; i < repeats; ++i) {
        const auto start = Clock::now();
        const double value = f();
        const auto end = Clock::now();
        sink = sink + value;
        times.push_back(std::chrono::duration<double>(end - start).count());
    }
    std::sort(times.begin(), times.end());
    return times[times.size()/2];
}
}
int main(int argc, char** argv) {
    if (argc < 2 || argc > 5) {
        std::cerr << "usage: kumo-automatic-backend-benchmark FILE.fchk[.gz] [points=512] [block=128] [repeats=5]\n";
        return 2;
    }
    try {
        // Loading and orbital eligibility validation are excluded from timing.
        const auto molecule = Molecule::from_fchk(argv[1]);
        if (!molecule.electronic_structure_info().ao_mo_eligible())
            throw std::runtime_error("AO_MO is not eligible for this FCHK");
        const int count = argc > 2 ? std::stoi(argv[2]) : 512;
        const int block = argc > 3 ? std::stoi(argv[3]) : 128;
        const int repeats = argc > 4 ? std::stoi(argv[4]) : 5;
        if (count <= 0 || block <= 0 || repeats <= 0)
            throw std::invalid_argument("points, block and repeats must be positive");
        Grid grid;
        grid.origin = Eigen::Vector3d(-1.8, -1.2, -0.9);
        grid.steps = Eigen::Matrix3d::Identity() * 0.12;
        grid.counts = Eigen::Vector3i(count, 1, 1);
        const Evaluator ao(molecule, EvaluationOptions{1, static_cast<std::size_t>(block), Backend::AO_P});
        const Evaluator mo(molecule, EvaluationOptions{1, static_cast<std::size_t>(block), Backend::AO_MO});
        const Evaluator automatic(molecule, EvaluationOptions{1, static_cast<std::size_t>(block), Backend::AUTO});
        const auto& info = molecule.electronic_structure_info();
        std::cout << std::setprecision(10);
        for (auto workload : {Workload::Density, Workload::ELF, Workload::Hessian}) {
            const auto req = request_for(workload);
            const auto plan = detail::EvaluationPlan::from(req);
            const Backend selected = detail::choose_backend(Backend::AUTO,
                info.ao_mo_eligible(), molecule.basis_size(),
                info.occupied_alpha, info.occupied_beta, plan);
            const auto result_ao = ao.evaluate(grid, req);
            const auto result_mo = mo.evaluate(grid, req);
            const auto result_auto = automatic.evaluate(grid, req);
            const auto comparison = compare(result_ao, result_mo, req);
            const auto selected_comparison = compare(
                result_auto, selected == Backend::AO_MO ? result_mo : result_ao, req);
            if (selected_comparison.max_abs > 1e-12)
                throw std::runtime_error("AUTO result differs from selected explicit backend");
            const auto ao_run = [&] { return checksum(ao.evaluate(grid, req), req); };
            const auto mo_run = [&] { return checksum(mo.evaluate(grid, req), req); };
            const auto auto_run = [&] { return checksum(automatic.evaluate(grid, req), req); };
            (void)ao_run(); (void)mo_run(); (void)auto_run();
            // Reverse order in the second pass to reduce order-dependent bias.
            const double ta1 = median_seconds(ao_run, repeats);
            const double tm1 = median_seconds(mo_run, repeats);
            const double tm2 = median_seconds(mo_run, repeats);
            const double ta2 = median_seconds(ao_run, repeats);
            const double tu1 = median_seconds(auto_run, repeats);
            const double tu2 = median_seconds(auto_run, repeats);
            const double tu = (tu1 + tu2) / 2.0;
            const double ta = (ta1 + ta2) / 2.0;
            const double tm = (tm1 + tm2) / 2.0;
            std::cout << "workload=" << name(workload)
                      << " points=" << count << " block=" << block
                      << " repeats=" << repeats
                      << " selected=" << backend_name(selected)
                      << " ao_p_seconds=" << ta << " ao_mo_seconds=" << tm
                      << " auto_seconds=" << tu
                      << " auto_over_best=" << tu / std::min(ta, tm)
                      << " selection_regret=" << (selected == Backend::AO_MO ? tm : ta) / std::min(ta, tm)
                      << " speedup=" << ta/tm
                      << " max_abs=" << comparison.max_abs
                      << " max_scaled=" << comparison.max_scaled
                      << " checksum_ao=" << comparison.checksum_ao
                      << " checksum_mo=" << comparison.checksum_mo << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "automatic backend benchmark: " << e.what() << '\n';
        return 1;
    }
}
