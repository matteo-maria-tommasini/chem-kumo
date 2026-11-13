#include "chem/kumo.hpp"
#include "omp_support.hpp"
#include "test_support.hpp"
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
            const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 137);
            for (Eigen::Index i = 0; i < points.cols(); ++i) {
                points.col(i) =
                    Eigen::Vector3d(0.11 + 0.073 * i, -1.3 + 0.041 * i, 0.27 - 0.037 * i);
            }
            kumo_test::set_threads(1);
            for (auto spin : {chem::kumo::Spin::Total, chem::kumo::Spin::Alpha,
                              chem::kumo::Spin::Beta, chem::kumo::Spin::Difference}) {
                const auto reference = mol.density_derivatives_blocked(points, spin, 32);
                for (int threads : {1, 2, 4}) {
                    kumo_test::set_threads(threads);
                    for (std::size_t block :
                         {std::size_t(1), std::size_t(16), std::size_t(32), std::size_t(64)}) {
                        const auto actual = mol.density_derivatives_blocked(points, spin, block);
                        const std::string label = std::string(argv[file]) + " spin " +
                                                  std::to_string(static_cast<int>(spin)) +
                                                  " threads " + std::to_string(threads) +
                                                  " block " + std::to_string(block);
                        if (!kumo_test::check(actual.size() == reference.size(),
                                              label + ": same number of points as the reference"))
                            continue;
                        std::size_t mismatches = 0;
                        for (std::size_t i = 0; i < actual.size(); ++i) {
                            const auto& a = actual[i];
                            const auto& b = reference[i];
                            if (!std::isfinite(a.value) || !a.gradient.allFinite() ||
                                !a.hessian.allFinite() || !std::isfinite(a.laplacian) ||
                                !std::isfinite(a.tau) || std::abs(a.value - b.value) > 1e-10 ||
                                (a.gradient - b.gradient).cwiseAbs().maxCoeff() > 1e-10 ||
                                (a.hessian - b.hessian).cwiseAbs().maxCoeff() > 1e-10 ||
                                std::abs(a.laplacian - b.laplacian) > 1e-10 ||
                                std::abs(a.tau - b.tau) > 1e-10)
                                ++mismatches;
                        }
                        kumo_test::check(mismatches == 0,
                                         label +
                                             ": value, gradient, Hessian, Laplacian and tau finite "
                                             "and within 1e-10 of the 1-thread, block-32 "
                                             "reference (" +
                                             std::to_string(mismatches) + " points differ)");
                    }
                }
            }
            if (kumo_test::failures == failures_before)
                std::cout << argv[file] << " parallel/serial derivatives OK\n";
        }
    });
}
