#include "chem/kumo.hpp"
#include "omp_support.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    return kumo_test::run([&] {
        for (int file = 1; file < argc; ++file) {
            const int failures_before = kumo_test::failures;
            const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
            const auto grid = chem::kumo::Grid::around(mol.atoms(), 2.0, 0.55);
            const auto components = {chem::kumo::Spin::Total, chem::kumo::Spin::Alpha,
                                     chem::kumo::Spin::Beta, chem::kumo::Spin::Difference};
            kumo_test::set_threads(1);
            const auto reference = mol.densities_blocked(grid, components, 64);
            // Single-thread results by block size: with the same blocks, other
            // thread counts must give identical values (blocks are independent).
            std::map<std::size_t, std::map<chem::kumo::Spin, chem::kumo::Field>> serial;
            for (int threads : {1, 2, 4}) {
                kumo_test::set_threads(threads);
                for (std::size_t block :
                     {std::size_t(1), std::size_t(32), std::size_t(64), std::size_t(128)}) {
                    const auto fields = mol.densities_blocked(grid, components, block);
                    if (threads == 1) serial.emplace(block, fields);
                    for (auto spin : components) {
                        const auto a = reference.at(spin).values();
                        const auto b = fields.at(spin).values();
                        const std::string label = std::string(argv[file]) + " spin " +
                                                  std::to_string(static_cast<int>(spin)) +
                                                  " threads " + std::to_string(threads) +
                                                  " block " + std::to_string(block);
                        if (!kumo_test::check(a.size() == b.size(),
                                              label + ": same number of grid points"))
                            continue;
                        std::size_t mismatches = 0;
                        double worst = 0.0;
                        for (std::size_t i = 0; i < a.size(); ++i) {
                            if (!std::isfinite(b[i]) || std::abs(a[i] - b[i]) > 1e-10) {
                                ++mismatches;
                                worst = std::max(worst, std::abs(a[i] - b[i]));
                            }
                        }
                        kumo_test::check(mismatches == 0,
                                         label +
                                             ": finite and within 1e-10 of the 1-thread, block-64 "
                                             "reference (" +
                                             std::to_string(mismatches) +
                                             " points differ, max_abs=" + kumo_test::num(worst) +
                                             ")");
                        const auto same = serial.at(block).at(spin).values();
                        kumo_test::check(std::equal(same.begin(), same.end(), b.begin(), b.end()),
                                         label +
                                             ": identical to 1 thread with the same block size");
                    }
                }
            }
            if (kumo_test::failures == failures_before)
                std::cout << argv[file] << " parallel/serial density comparison OK\n";
        }
    });
}
