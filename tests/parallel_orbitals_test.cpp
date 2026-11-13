#include "chem/kumo.hpp"
#include "omp_support.hpp"
#include "test_support.hpp"
#include <map>
#include <vector>
#include <Eigen/Core>
#include <algorithm>
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
            const std::string name = argv[file];
            const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
            const auto grid = chem::kumo::Grid::around(mol.atoms(), 2.0, 0.55);
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 137);
            for (Eigen::Index j = 0; j < points.cols(); ++j)
                points.col(j) =
                    Eigen::Vector3d(0.21 + 0.043 * j, -0.37 + 0.017 * j, 0.14 - 0.031 * j);
            kumo_test::set_threads(1);
            for (const auto spin : {chem::kumo::Spin::Alpha, chem::kumo::Spin::Beta}) {
                const auto n = mol.orbital_count(spin);
                if (!kumo_test::check(n != 0, name + ": molecular orbitals of spin " +
                                                  std::to_string(static_cast<int>(spin)) +
                                                  " are present"))
                    continue;
                const std::vector<std::size_t> indices = {
                    1, std::min<std::size_t>(3, n),
                    std::min<std::size_t>(static_cast<std::size_t>(mol.occupied_count(spin)) + 1,
                                          n)};
                const auto reference = mol.orbitals(grid, indices, spin);
                const auto derivative_reference =
                    mol.orbital_derivatives_blocked(points, indices[1], spin, 64);
                // Single-thread derivatives by block size: with the same blocks,
                // other thread counts must give identical values.
                std::map<std::size_t, std::vector<chem::kumo::OrbitalDerivatives>> serial;
                for (int threads : {1, 2, 4}) {
                    kumo_test::set_threads(threads);
                    const std::string label = name + " spin " +
                                              std::to_string(static_cast<int>(spin)) + " threads " +
                                              std::to_string(threads);
                    const auto fields = mol.orbitals(grid, indices, spin);
                    if (kumo_test::check(fields.size() == reference.size(),
                                         label + ": same number of MO fields as the reference")) {
                        for (std::size_t k = 0; k < fields.size(); ++k) {
                            const auto a = reference[k].values(), b = fields[k].values();
                            const std::string mo = label + " MO " + std::to_string(indices[k]);
                            if (!kumo_test::check(a.size() == b.size(),
                                                  mo + ": same number of grid points"))
                                continue;
                            // Same blocks as the reference, so identical values.
                            std::size_t mismatches = 0;
                            for (std::size_t i = 0; i < a.size(); ++i)
                                if (!std::isfinite(b[i]) || a[i] != b[i]) ++mismatches;
                            kumo_test::check(mismatches == 0,
                                             mo +
                                                 ": grid values finite and identical to the "
                                                 "1-thread reference (" +
                                                 std::to_string(mismatches) + " points differ)");
                        }
                    }
                    for (std::size_t block :
                         {std::size_t(1), std::size_t(32), std::size_t(64), std::size_t(128)}) {
                        const auto derivatives =
                            mol.orbital_derivatives_blocked(points, indices[1], spin, block);
                        if (threads == 1) serial.emplace(block, derivatives);
                        const std::string mo = label + " MO " + std::to_string(indices[1]) +
                                               " block " + std::to_string(block);
                        if (!kumo_test::check(derivatives.size() == derivative_reference.size(),
                                              mo + ": same number of derivative points"))
                            continue;
                        std::size_t mismatches = 0;
                        for (std::size_t j = 0; j < derivatives.size(); ++j) {
                            const auto& a = derivative_reference[j];
                            const auto& b = derivatives[j];
                            if (!std::isfinite(b.value) || !b.gradient.allFinite() ||
                                !b.hessian.allFinite() || std::abs(a.value - b.value) > 1e-10 ||
                                (a.gradient - b.gradient).cwiseAbs().maxCoeff() > 1e-10 ||
                                (a.hessian - b.hessian).cwiseAbs().maxCoeff() > 1e-10)
                                ++mismatches;
                        }
                        kumo_test::check(mismatches == 0,
                                         mo +
                                             ": value, gradient and Hessian finite and within "
                                             "1e-10 of the 1-thread, block-64 reference (" +
                                             std::to_string(mismatches) + " points differ)");
                        const auto& same = serial.at(block);
                        std::size_t different = 0;
                        for (std::size_t j = 0; j < derivatives.size(); ++j)
                            if (derivatives[j].value != same[j].value ||
                                derivatives[j].gradient != same[j].gradient ||
                                derivatives[j].hessian != same[j].hessian)
                                ++different;
                        kumo_test::check(different == 0,
                                         mo + ": identical to 1 thread with the same block size (" +
                                             std::to_string(different) + " points differ)");
                    }
                }
            }
            kumo_test::check_throws<std::invalid_argument>(
                [&] {
                    (void)mol.orbital_derivatives_blocked(points, 1, chem::kumo::Spin::Alpha, 0);
                },
                name + ": orbital_derivatives_blocked rejects block size 0");
            if (kumo_test::failures == failures_before)
                std::cout << argv[file] << " parallel orbital grid/derivatives OK\n";
        }
    });
}
