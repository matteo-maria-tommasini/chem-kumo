#include "chem/kumo.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
namespace {
using Clock = std::chrono::steady_clock;
volatile double sink = 0.0;
struct Measurement { double seconds; double checksum; };

template <class F> Measurement measure(F&& function, int repeats) {
    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(repeats));
    double checksum = 0.0;
    for (int i = 0; i < repeats; ++i) {
        const auto start = Clock::now();
        checksum = function();
        const auto stop = Clock::now();
        sink = sink + checksum;
        times.push_back(std::chrono::duration<double>(stop - start).count());
    }
    std::sort(times.begin(), times.end());
    return {times[times.size() / 2], checksum};
}

constexpr std::array components{DensityComponent::Total, DensityComponent::Alpha,
                                DensityComponent::Beta, DensityComponent::Spin};
Request make_request(bool hessian) {
    Request request;
    for (const auto component : components) {
        const DensityDescriptor density{component};
        request.add(density.gradient()).add(density.tau());
        if (hessian) request.add(density.hessian());
    }
    return request;
}

double checksum_common(const Results& result) {
    double checksum = 0.0;
    for (const auto component : components) {
        const DensityDescriptor density{component};
        const auto& gradient = result.get(density.gradient());
        const auto& tau = result.get(density.tau());
        for (std::size_t p = 0; p < tau.points; ++p)
            checksum += gradient.component(0)[p] + gradient.component(1)[p]
                      + gradient.component(2)[p] + tau.component()[p];
    }
    return checksum;
}
void validate(const Results& low, const Results& full) {
    for (const auto component : components) {
        const DensityDescriptor density{component};
        for (const auto property : {density.gradient(), density.tau()}) {
            const auto& a = low.get(property);
            const auto& b = full.get(property);
            if (a.points != b.points || a.components.size() != b.components.size())
                throw std::runtime_error("different result dimensions");
            for (std::size_t c = 0; c < a.components.size(); ++c)
                for (std::size_t p = 0; p < a.points; ++p) {
                    const double x = a.component(c)[p], y = b.component(c)[p];
                    if (!std::isfinite(x) || !std::isfinite(y) ||
                        std::abs(x-y) > 5e-11 * std::max({1.0, std::abs(x), std::abs(y)}))
                        throw std::runtime_error("first/second-order numerical mismatch");
                }
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 5) {
        std::cerr << "usage: kumo-selective-ao-benchmark FILE.fchk[.gz] [points=512] [block=128] [repeats=5]\n";
        return 2;
    }
    try {
        const auto molecule = Molecule::from_fchk(argv[1]);
        const int count = argc > 2 ? std::stoi(argv[2]) : 512;
        const int block = argc > 3 ? std::stoi(argv[3]) : 128;
        const int repeats = argc > 4 ? std::stoi(argv[4]) : 5;
        if (count < 1 || block < 1 || repeats < 1)
            throw std::invalid_argument("positive points, block, and repeats required");
        Grid grid;
        grid.origin = Eigen::Vector3d(-1.8, -1.2, -0.9);
        grid.steps = Eigen::Matrix3d::Identity() * 0.12;
        grid.counts = Eigen::Vector3i(count, 1, 1);
        const Evaluator evaluator(molecule, EvaluationOptions{1, static_cast<std::size_t>(block)});
        const Request first = make_request(false);
        const Request second = make_request(true);
        const auto low = evaluator.evaluate(grid, first);
        const auto full = evaluator.evaluate(grid, second);
        validate(low, full);
        const auto first_order = [&] { return checksum_common(evaluator.evaluate(grid, first)); };
        const auto second_order = [&] { return checksum_common(evaluator.evaluate(grid, second)); };
        // Warm both paths; alternate measurement order across two passes.
        (void)first_order(); (void)second_order();
        const auto a = measure(first_order, repeats);
        const auto b = measure(second_order, repeats);
        const auto c = measure(second_order, repeats);
        const auto d = measure(first_order, repeats);
        const double first_seconds = 0.5 * (a.seconds + d.seconds);
        const double second_seconds = 0.5 * (b.seconds + c.seconds);
        std::cout << "points=" << count << " block=" << block << " repeats=" << repeats
                  << " first_order_seconds=" << first_seconds
                  << " second_order_seconds=" << second_seconds
                  << " ratio=" << second_seconds / first_seconds
                  << " checksum=" << a.checksum << '\n';
        // The second-order request additionally materializes Hessians. This is
        // an end-to-end workload comparison, not an isolated AO-kernel speedup.
    } catch (const std::exception& e) {
        std::cerr << "selective AO benchmark: " << e.what() << '\n';
        return 1;
    }
}
