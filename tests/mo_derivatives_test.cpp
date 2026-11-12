#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    constexpr double h = 1.e-4;
    return kumo_test::run([&] {
        for (int file = 1; file < argc; ++file) {
            const std::string name = argv[file];
            const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
            double max_gradient = 0.0, max_hessian = 0.0;
            for (const auto channel : {chem::kumo::Spin::Alpha, chem::kumo::Spin::Beta}) {
                const auto count = mol.orbital_count(channel);
                if (!kumo_test::check(count != 0, name + ": molecular orbitals of spin " +
                                                      std::to_string(static_cast<int>(channel)) +
                                                      " are present"))
                    continue;
                const std::size_t indices[] = {
                    1, std::min<std::size_t>(3, count),
                    std::min<std::size_t>(static_cast<std::size_t>(mol.occupied_count(channel)) + 1,
                                          count)};
                for (const auto index : indices) {
                    for (int t = 0; t < 8; ++t) {
                        const Eigen::Vector3d r(0.41 + 0.33 * t, -0.37 + 0.17 * t, 0.23 - 0.13 * t);
                        const auto exact = mol.orbital_derivatives(r, index, channel);
                        // The analytic Hessian is assembled symmetrically; 1e-13
                        // allows rounding only.
                        kumo_test::check(
                            std::isfinite(exact.value) && exact.gradient.allFinite() &&
                                exact.hessian.allFinite() &&
                                !((exact.hessian - exact.hessian.transpose())
                                      .cwiseAbs()
                                      .maxCoeff() > 1.e-13),
                            name + " MO " + std::to_string(index) + " spin " +
                                std::to_string(static_cast<int>(channel)) + " point " +
                                std::to_string(t) +
                                ": analytic derivatives finite, Hessian symmetric within 1e-13");
                        for (int a = 0; a < 3; ++a) {
                            Eigen::Vector3d rp = r, rm = r;
                            rp[a] += h;
                            rm[a] -= h;
                            const auto plus = mol.orbital_derivatives(rp, index, channel);
                            const auto minus = mol.orbital_derivatives(rm, index, channel);
                            const double fd1 = (plus.value - minus.value) / (2 * h);
                            max_gradient =
                                std::max(max_gradient, std::abs(fd1 - exact.gradient[a]));
                            for (int b = 0; b < 3; ++b) {
                                const double fd2 = (plus.gradient[b] - minus.gradient[b]) / (2 * h);
                                max_hessian =
                                    std::max(max_hessian, std::abs(fd2 - exact.hessian(b, a)));
                            }
                        }
                    }
                }
            }
            std::cout << argv[file] << " grad=" << max_gradient << " hessian=" << max_hessian
                      << '\n';
            // Central differences with h = 1e-4 have truncation error
            // h^2/6 f''' (plus rounding ~eps |f| / h); third derivatives of
            // tight, high-exponent basis functions are large, hence the loose
            // bounds (1e-3 gradient, 1e-2 Hessian). They still catch wrong
            // derivative formulas, whose errors are much larger.
            kumo_test::check(!(max_gradient > 1.e-3),
                             name +
                                 ": MO gradient matches central differences within 1e-3 "
                                 "(max_abs=" +
                                 kumo_test::num(max_gradient) + ")");
            kumo_test::check(!(max_hessian > 1.e-2),
                             name +
                                 ": MO Hessian matches central differences within 1e-2 "
                                 "(max_abs=" +
                                 kumo_test::num(max_hessian) + ")");
            kumo_test::check_throws<std::out_of_range>(
                [&] { (void)mol.orbital_derivatives(Eigen::Vector3d::Zero(), 0); },
                name + ": orbital_derivatives rejects MO index 0");
            kumo_test::check_throws<std::invalid_argument>(
                [&] {
                    (void)mol.orbital_derivatives(Eigen::Vector3d::Zero(), 1,
                                                  chem::kumo::Spin::Total);
                },
                name + ": orbital_derivatives rejects the Total spin channel");
        }
    });
}
