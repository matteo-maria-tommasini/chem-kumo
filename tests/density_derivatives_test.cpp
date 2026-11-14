#include "chem/kumo.hpp"
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
std::string context; // file and point of the current checks

void require(bool condition, const std::string& message) {
    kumo_test::check(condition, context + message);
}

// Both finite and |a - b| <= absolute + relative * max(|a|, |b|). (Not
// kumo_test::check_close: the tolerance is a sum, not a maximum.)
void close(double a, double b, double absolute, double relative, const char* message) {
    const double tolerance = absolute + relative * std::max(std::abs(a), std::abs(b));
    require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tolerance,
            std::string(message) + " (got " + kumo_test::num(a, 17) + " and " +
                kumo_test::num(b, 17) + ", tolerance " + kumo_test::num(tolerance) + ")");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    constexpr double step = 1.e-4;
    return kumo_test::run([&] {
        for (int file = 1; file < argc; ++file) {
            const auto molecule = chem::kumo::Molecule::from_fchk(argv[file]);
            double max_gradient = 0.0, max_hessian = 0.0;
            for (int t = 0; t < 9; ++t) {
                const Eigen::Vector3d r(0.41 + 0.29 * t, -0.37 + 0.13 * t, 0.23 - 0.17 * t);
                context = std::string(argv[file]) + " point " + std::to_string(t) + ": ";
                const auto total = kumo_test::Kernels::density_derivatives(molecule, r);
                const auto alpha = kumo_test::Kernels::density_derivatives(
                    molecule, r, chem::kumo::DensityComponent::Alpha);
                const auto beta = kumo_test::Kernels::density_derivatives(
                    molecule, r, chem::kumo::DensityComponent::Beta);
                const auto spin = kumo_test::Kernels::density_derivatives(
                    molecule, r, chem::kumo::DensityComponent::Difference);
                for (const auto* f : {&total, &alpha, &beta, &spin}) {
                    require(std::isfinite(f->value) && f->gradient.allFinite() &&
                                f->hessian.allFinite() && std::isfinite(f->tau),
                            "density derivatives are finite");
                    close(f->laplacian, f->hessian.trace(), 1.e-12, 1.e-12,
                          "Laplacian equals the Hessian trace");
                    require((f->hessian - f->hessian.transpose()).cwiseAbs().maxCoeff() < 1.e-12,
                            "Hessian symmetric within 1e-12");
                }
                for (int a = 0; a < 3; ++a) {
                    Eigen::Vector3d rp = r, rm = r;
                    rp[a] += step;
                    rm[a] -= step;
                    for (auto component :
                         {chem::kumo::DensityComponent::Total, chem::kumo::DensityComponent::Alpha,
                          chem::kumo::DensityComponent::Beta,
                          chem::kumo::DensityComponent::Difference}) {
                        const auto exact =
                            kumo_test::Kernels::density_derivatives(molecule, r, component);
                        const auto plus =
                            kumo_test::Kernels::density_derivatives(molecule, rp, component);
                        const auto minus =
                            kumo_test::Kernels::density_derivatives(molecule, rm, component);
                        const double fd = (plus.value - minus.value) / (2 * step);
                        max_gradient = std::max(max_gradient, std::abs(fd - exact.gradient[a]));
                        close(fd, exact.gradient[a], 2.e-4, 2.e-5,
                              "density gradient matches central differences");
                        for (int b = 0; b < 3; ++b) {
                            const double second =
                                (plus.gradient[b] - minus.gradient[b]) / (2 * step);
                            max_hessian =
                                std::max(max_hessian, std::abs(second - exact.hessian(b, a)));
                            close(second, exact.hessian(b, a), 3.e-3, 3.e-5,
                                  "density Hessian matches central differences");
                        }
                    }
                }
                // All fields, including tau, must obey the same spin decomposition.
                close(total.value, alpha.value + beta.value, 1.e-10, 1.e-12,
                      "total density = alpha + beta");
                close(spin.value, alpha.value - beta.value, 1.e-10, 1.e-12,
                      "spin density = alpha - beta");
                close(total.laplacian, alpha.laplacian + beta.laplacian, 1.e-9, 1.e-12,
                      "total Laplacian = alpha + beta");
                close(spin.laplacian, alpha.laplacian - beta.laplacian, 1.e-9, 1.e-12,
                      "spin Laplacian = alpha - beta");
                close(total.tau, alpha.tau + beta.tau, 1.e-10, 1.e-12, "total tau = alpha + beta");
                close(spin.tau, alpha.tau - beta.tau, 1.e-10, 1.e-12, "spin tau = alpha - beta");
                require((total.gradient - alpha.gradient - beta.gradient).norm() < 1.e-9,
                        "total gradient = alpha + beta within 1e-9");
                require((spin.gradient - alpha.gradient + beta.gradient).norm() < 1.e-9,
                        "spin gradient = alpha - beta within 1e-9");
                require((total.hessian - alpha.hessian - beta.hessian).norm() < 1.e-8,
                        "total Hessian = alpha + beta within 1e-8");
                require((spin.hessian - alpha.hessian + beta.hessian).norm() < 1.e-8,
                        "spin Hessian = alpha - beta within 1e-8");
                // Independent pointwise density implementation, not derivative code.
                chem::kumo::Grid one;
                one.origin = r;
                one.counts = Eigen::Vector3i::Ones();
                close(total.value, kumo_test::Kernels::density_pointwise(molecule, one)(0, 0, 0),
                      1.e-10, 1.e-12, "density value matches density_pointwise()");
            }
            kumo_test::check_throws<std::invalid_argument>(
                [&] {
                    (void)kumo_test::Kernels::density_derivatives(molecule,
                                                                  Eigen::Vector3d::Constant(NAN));
                },
                std::string(argv[file]) + ": density_derivatives rejects a nonfinite point");
            std::cout << argv[file] << " grad=" << max_gradient << " hessian=" << max_hessian
                      << '\n';
        }
    });
}
