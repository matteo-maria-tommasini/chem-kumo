#include "chem/kumo.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
using Clock = std::chrono::steady_clock;

chem::kumo::Grid make_grid(const std::vector<chem::kumo::Atom>& atoms,
                           std::size_t requested, double margin) {
    if (atoms.empty()) throw std::invalid_argument("empty molecule");
    const auto side = std::ceil(std::cbrt(static_cast<double>(requested)));
    if (side < 2.0 || side > 1000.0)
        throw std::invalid_argument("invalid grid size");
    const int n = static_cast<int>(side);
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(n, n, n);
    Eigen::Vector3d lo = atoms.front().position;
    Eigen::Vector3d hi = lo;
    for (const auto& atom : atoms) {
        lo = lo.cwiseMin(atom.position);
        hi = hi.cwiseMax(atom.position);
    }
    lo.array() -= margin;
    hi.array() += margin;
    grid.origin = lo;
    grid.steps.setZero();
    for (int k = 0; k < 3; ++k)
        grid.steps(k, k) = (hi[k] - lo[k]) / (n - 1);
    return grid;
}

struct Result {
    std::vector<double> total;
    std::vector<double> spin;
};

Result evaluate(const chem::kumo::Molecule& mol, const chem::kumo::Grid& grid,
                std::size_t block_size) {
    auto fields = mol.densities_blocked(grid,
        {chem::kumo::Spin::Total, chem::kumo::Spin::Difference}, block_size);
    const auto a = fields.at(chem::kumo::Spin::Total).values();
    const auto b = fields.at(chem::kumo::Spin::Difference).values();
    return {std::vector<double>(a.begin(), a.end()),
            std::vector<double>(b.begin(), b.end())};
}

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::runtime_error("grid size mismatch");
    double diff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i]))
            throw std::runtime_error("nonfinite density sample");
        diff = std::max(diff, std::abs(a[i] - b[i]));
    }
    return diff;
}

void set_threads(int threads) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#else
    if (threads != 1) throw std::runtime_error("OpenMP not enabled");
#endif
}

template <typename F>
double measure(int repeats, F&& work) {
    work(); // warm-up, excluded
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < repeats; ++i) {
        const auto start = Clock::now();
        work();
        const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
        best = std::min(best, seconds);
    }
    return best;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 6) {
        std::cerr << "Usage: kumo-parallel-density-benchmark FILE.fchk[.gz] "
                     "[points=10000] [repeats=3] [block_size=128] [margin_bohr=2.5]\n";
        return 2;
    }
    try {
        const std::size_t requested = argc > 2 ? std::stoull(argv[2]) : 10000;
        const int repeats = argc > 3 ? std::stoi(argv[3]) : 3;
        const std::size_t block_size = argc > 4 ? std::stoull(argv[4]) : 128;
        const double margin = argc > 5 ? std::stod(argv[5]) : 2.5;
        if (requested < 8 || requested > 1000000 || repeats < 1 || repeats > 20 ||
            block_size < 1 || block_size > 65536 || !std::isfinite(margin) || margin <= 0)
            throw std::invalid_argument("invalid benchmark arguments");

        // Parsing and construction are excluded from every timed section.
        const auto mol = chem::kumo::Molecule::from_fchk(argv[1]);
        const auto grid = make_grid(mol.atoms(), requested, margin);
        const std::size_t points = static_cast<std::size_t>(grid.counts[0]) *
            static_cast<std::size_t>(grid.counts[1]) * static_cast<std::size_t>(grid.counts[2]);
        std::cout << std::setprecision(12)
                  << "file=" << argv[1] << "\n"
                  << "basis=" << mol.basis_size() << " atoms=" << mol.atoms().size()
                  << " points=" << points << " repeats=" << repeats
                  << " block_size=" << block_size << " margin_bohr=" << margin << "\n"
#ifdef _OPENMP
                  << "openmp=enabled\n";
#else
                  << "openmp=disabled\n";
#endif
        set_threads(1);
        // Independent, single-thread reference computed before the sweep.
        const auto reference = evaluate(mol, grid, block_size);
        Result latest;
        const double baseline = measure(repeats, [&] {
            latest = evaluate(mol, grid, block_size);
        });
        const double reference_error = std::max(max_difference(reference.total, latest.total),
                                                 max_difference(reference.spin, latest.spin));
        if (reference_error > 1e-10)
            throw std::runtime_error("single-thread reference mismatch");
        const double checksum = reference.total.at(points / 2) + reference.spin.at(points / 2);
        std::cout << "threads=1 time_s=" << baseline
                  << " speedup=1 efficiency=1 max_abs=" << reference_error << "\n";

#ifdef _OPENMP
        for (const int threads : {2, 4, 8, 16}) {
            set_threads(threads);
            const double seconds = measure(repeats, [&] {
                latest = evaluate(mol, grid, block_size);
            });
            const double error = std::max(max_difference(reference.total, latest.total),
                                          max_difference(reference.spin, latest.spin));
            if (error > 1e-10)
                throw std::runtime_error("parallel density differs from serial reference");
            const double speedup = baseline / seconds;
            std::cout << "threads=" << threads << " time_s=" << seconds
                      << " speedup=" << speedup
                      << " efficiency=" << speedup / threads
                      << " max_abs=" << error << "\n";
        }
#endif
        std::cout << "checksum=" << checksum << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
