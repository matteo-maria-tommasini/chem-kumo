#include "chem/kumo.hpp"
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
// The Evaluator's OrbitalGradient and OrbitalHessian, at the points of a small
// grid and at the same points given as a matrix, against the pointwise
// reference: equal to rounding (1e-12 of the largest magnitude).
void check_evaluator(const chem::kumo::Molecule& mol, const std::string& name,
                     chem::kumo::SpinChannel channel, std::size_t index) {
    using namespace chem::kumo;
    const auto orbital =
        channel == SpinChannel::Alpha ? orbital::alpha(index) : orbital::beta(index);
    const Grid grid =
        Grid::plane(mol.atoms().front().position, mol.atoms().back().position,
                    mol.atoms().front().position + Eigen::Vector3d(0.3, 0.7, 0.2), 0.8, 0.45);
    const auto n = static_cast<Eigen::Index>(static_cast<std::size_t>(grid.counts.prod()));
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, n);
    for (Eigen::Index q = 0; q < n; ++q)
        points.col(q) =
            grid.origin +
            grid.steps * Eigen::Vector3d(static_cast<double>(q / grid.counts[2] / grid.counts[1]),
                                         static_cast<double>(q / grid.counts[2] % grid.counts[1]),
                                         static_cast<double>(q % grid.counts[2]));
    const auto request = Request().add(orbital.gradient()).add(orbital.hessian());
    const Evaluator evaluator(mol);
    const auto on_grid = evaluator.evaluate(grid, request);
    const auto at_points = evaluator.evaluate(points, request);
    constexpr int column[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    double worst = 0.0;
    for (const auto* results : {&on_grid, &at_points})
        for (Eigen::Index q = 0; q < n; ++q) {
            const auto p = static_cast<std::size_t>(q);
            const auto exact =
                kumo_test::Kernels::orbital_derivatives(mol, points.col(q), index, channel);
            const double scale = std::max(
                {1.0, exact.gradient.cwiseAbs().maxCoeff(), exact.hessian.cwiseAbs().maxCoeff()});
            for (std::size_t a = 0; a < 3; ++a) {
                worst = std::max(worst, std::abs(results->get(orbital.gradient()).component(a)[p] -
                                                 exact.gradient[static_cast<Eigen::Index>(a)]) /
                                            scale);
                for (int b = 0; b < 3; ++b)
                    worst = std::max(
                        worst, std::abs(results->get(orbital.hessian())
                                            .component(static_cast<std::size_t>(column[a][b]))[p] -
                                        exact.hessian(static_cast<Eigen::Index>(a), b)) /
                                   scale);
            }
        }
    kumo_test::check(worst <= 1e-12, name + " MO " + std::to_string(index) +
                                         ": Evaluator orbital gradient and Hessian, on a grid "
                                         "and at points, equal the pointwise reference to 1e-12 "
                                         "(worst " +
                                         kumo_test::num(worst) + ")");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    constexpr double h = 1.e-4;
    return kumo_test::run([&] {
        for (int file = 1; file < argc; ++file) {
            const std::string name = argv[file];
            const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
            double max_gradient = 0.0, max_hessian = 0.0;
            for (const auto channel :
                 {chem::kumo::SpinChannel::Alpha, chem::kumo::SpinChannel::Beta}) {
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
                    check_evaluator(mol, name, channel, index);
                    for (int t = 0; t < 8; ++t) {
                        const Eigen::Vector3d r(0.41 + 0.33 * t, -0.37 + 0.17 * t, 0.23 - 0.13 * t);
                        const auto exact =
                            kumo_test::Kernels::orbital_derivatives(mol, r, index, channel);
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
                            const auto plus =
                                kumo_test::Kernels::orbital_derivatives(mol, rp, index, channel);
                            const auto minus =
                                kumo_test::Kernels::orbital_derivatives(mol, rm, index, channel);
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
                [&] {
                    (void)kumo_test::Kernels::orbital_derivatives(mol, Eigen::Vector3d::Zero(), 0);
                },
                name + ": orbital_derivatives rejects MO index 0");
            kumo_test::check_throws<std::invalid_argument>(
                [&] {
                    (void)kumo_test::Kernels::orbital_derivatives(
                        mol, Eigen::Vector3d::Zero(), 1, static_cast<chem::kumo::SpinChannel>(7));
                },
                name + ": orbital_derivatives rejects an invalid spin channel");
        }
    });
}
