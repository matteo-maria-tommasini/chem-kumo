#include "chem/kumo/density.hpp"
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;

namespace {
constexpr double algorithm_tol = 2e-10;

std::string where(const char* file, DensityComponent spin, std::size_t block, std::size_t point) {
    return std::string(" file=") + file + " spin=" + std::to_string(static_cast<int>(spin)) +
           " block=" + std::to_string(block) + " point=" + std::to_string(point);
}

// AO-MO result vs the AO-P contraction of the reconstructed MO density
// matrix: both finite and within algorithm_tol * max(1, |actual|, |reference|).
void check(double actual, double reference, const std::string& label, const char* file,
           DensityComponent spin, std::size_t block, std::size_t point) {
    const double scale = std::max({1.0, std::abs(actual), std::abs(reference)});
    const double error = std::abs(actual - reference);
    if (std::isfinite(actual) && std::isfinite(reference) && !(error > algorithm_tol * scale))
        return;
    kumo_test::check_close(actual, reference, algorithm_tol, algorithm_tol,
                           "AO-MO matches AO-P(P_MO) for " + label +
                               where(file, spin, block, point));
}

void compare(const detail::DensityDerivatives& x, const detail::DensityDerivatives& y,
             const char* file, DensityComponent spin, std::size_t block, std::size_t point,
             int order) {
    const auto test = [&](double a, double b, const std::string& name) {
        check(a, b, name, file, spin, block, point);
    };
    test(x.value, y.value, "density");
    for (int a = 0; a < 3; ++a)
        test(x.gradient[a], y.gradient[a], std::string("gradient_") + "xyz"[a]);
    test(x.tau, y.tau, "tau");
    if (order == 1) {
        test(x.laplacian, 0.0, "first-order laplacian");
        return;
    }
    test(x.laplacian, y.laplacian, "laplacian");
    constexpr int pairs[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
    constexpr const char* names[6] = {"xx", "yy", "zz", "xy", "xz", "yz"};
    for (int k = 0; k < 6; ++k)
        test(x.hessian(pairs[k][0], pairs[k][1]), y.hessian(pairs[k][0], pairs[k][1]),
             std::string("hessian_") + names[k]);
    test(x.hessian.trace(), x.laplacian, "MO Hessian trace");
    test(y.hessian.trace(), y.laplacian, "AO-P Hessian trace");
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_mo(x.hessian);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_ao(y.hessian);
    if (!kumo_test::check(eig_mo.info() == Eigen::Success && eig_ao.info() == Eigen::Success,
                          "Hessian eigensolver succeeds" + where(file, spin, block, point)))
        return;
    const double bound =
        (x.hessian - y.hessian).norm(); // Weyl: |delta lambda_i| <= ||delta H||_2 <= ||delta H||_F
    for (int k = 0; k < 3; ++k) {
        const double diff = std::abs(eig_mo.eigenvalues()[k] - eig_ao.eigenvalues()[k]);
        test(eig_mo.eigenvalues()[k], eig_ao.eigenvalues()[k],
             "hessian_eigenvalue_" + std::to_string(k));
        const double weyl =
            bound + 1e-10 * std::max({1.0, eig_mo.eigenvalues().cwiseAbs().maxCoeff(),
                                      eig_ao.eigenvalues().cwiseAbs().maxCoeff()});
        kumo_test::check(!(diff > weyl), "Hessian eigenvalue " + std::to_string(k) +
                                             " difference " + kumo_test::num(diff) +
                                             " within the Weyl bound " + kumo_test::num(weyl) +
                                             where(file, spin, block, point));
    }
    const double lambda2_mo = eig_mo.eigenvalues()[1];
    const double lambda2_ao = eig_ao.eigenvalues()[1];
    // Only compare signs when both are well separated from the error bound.
    kumo_test::check(!(std::min(std::abs(lambda2_mo), std::abs(lambda2_ao)) > bound + 1e-10 &&
                       std::signbit(lambda2_mo) != std::signbit(lambda2_ao)),
                     "stable lambda2 has the same sign in AO-MO (" + kumo_test::num(lambda2_mo) +
                         ") and AO-P (" + kumo_test::num(lambda2_ao) + ")" +
                         where(file, spin, block, point));
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        const std::vector<DensityComponent> spins{DensityComponent::Total, DensityComponent::Alpha,
                                                  DensityComponent::Beta,
                                                  DensityComponent::Difference};
        Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 19);
        for (int j = 0; j < 19; ++j)
            points.col(j) = Eigen::Vector3d(0.17 * j - 1.2, 0.09 * j - 0.8, 0.13 * j + 0.2);
        for (int arg = 1; arg < argc; ++arg) {
            const auto m = Molecule::from_fchk(argv[arg]);
            if (!kumo_test::check(m.electronic_structure_info().ao_mo_eligible(),
                                  std::string(argv[arg]) + ": the test molecule is AO-MO eligible"))
                continue;
            const auto matrices = detail::MoleculeTestAccess::reconstructed_matrices(m);
            for (int order : {1, 2})
                for (std::size_t block : {std::size_t{1}, std::size_t{7}, std::size_t{32}}) {
                    const auto mo = detail::MoleculeTestAccess::density_derivatives_mo(
                        m, points, spins, block, order);
                    for (DensityComponent spin : spins)
                        for (std::size_t j = 0; j < static_cast<std::size_t>(points.cols()); ++j) {
                            const auto ref = detail::MoleculeTestAccess::ao_p_reference(
                                m, points.col(static_cast<Eigen::Index>(j)), matrices.at(spin));
                            compare(mo.at(spin).at(j), ref, argv[arg], spin, block, j, order);
                        }
                }
        }
    });
}
