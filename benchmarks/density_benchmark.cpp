#include "chem/kumo.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: kumo-benchmark input.fchk [spacing_bohr=0.3] [repeats=3]\n";
        return 2;
    }
    using namespace chem::kumo;
    const double spacing = argc > 2 ? std::stod(argv[2]) : 0.3;
    const int repeats = argc > 3 ? std::stoi(argv[3]) : 3;
    if (repeats < 1 || repeats > 100) return 2;
    auto molecule = Molecule::from_fchk(argv[1]);
    auto grid = Grid::around(molecule.atoms(), 3.0, spacing);
    auto measure = [&](auto&& fn) {
        double best = 1e100;
        for (int k = 0; k < repeats; ++k) {
            const auto t0 = std::chrono::steady_clock::now();
            fn();
            const auto t1 = std::chrono::steady_clock::now();
            best = std::min(best, std::chrono::duration<double>(t1-t0).count());
        }
        return best;
    };
    volatile double sink = 0.0;
    const double separate = measure([&] {
        for (auto s : {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference}) {
            auto f = molecule.density(grid, s);
            sink = sink + f.values()[0];
        }
    });
    const double combined = measure([&] {
        auto f = molecule.spin_densities(grid);
        sink = sink + f.total.values()[0];
    });
    const double selective = measure([&] {
        auto f = molecule.densities(grid, {Spin::Total, Spin::Difference});
        sink = sink + f.at(Spin::Total).values()[0];
    });
    std::cout << "basis=" << molecule.basis_size()
              << " grid=" << static_cast<long long>(grid.counts[0])*grid.counts[1]*grid.counts[2]
              << " repeats=" << repeats << "\n"
              << "four_separate_s=" << separate << "\n"
              << "four_combined_s=" << combined << "\n"
              << "two_selected_s=" << selective << "\n"
              << "combined_speedup=" << separate/combined << "\n"
              << "sink=" << sink << "\n";
}
