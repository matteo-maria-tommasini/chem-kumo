#include "chem/kumo.hpp"
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
// The SCF density matrix and the MO coefficients are read from the FCHK
// file, which prints reals with 9 significant digits (E16.8), so densities
// built from the MOs and from the stored density matrix agree only to about
// that precision: 2e-5 absolute plus 2e-5 relative leaves a wide margin.
// (Not kumo_test::check_close: the tolerance is a sum, not a maximum.)
void check_close(double a, double b, const std::string& label) {
    const double tolerance = 2.e-5 + 2.e-5 * std::max(std::abs(a), std::abs(b));
    kumo_test::check(std::isfinite(a) && std::isfinite(b) && !(std::abs(a - b) > tolerance),
                     label + " (got " + kumo_test::num(a, 17) + " and " + kumo_test::num(b, 17) +
                         ", tolerance " + kumo_test::num(tolerance) + ")");
}

void compare(const chem::kumo::detail::DensityDerivatives& a,
             const chem::kumo::detail::DensityDerivatives& b, const std::string& where) {
    check_close(a.value, b.value, where + ": MO/DM densities agree");
    check_close(a.tau, b.tau, where + ": MO/DM tau agree");
    check_close(a.laplacian, b.laplacian, where + ": MO/DM Laplacians agree");
    for (int i = 0; i < 3; ++i) {
        check_close(a.gradient[i], b.gradient[i], where + ": MO/DM gradients agree");
        for (int j = 0; j < 3; ++j)
            check_close(a.hessian(i, j), b.hessian(i, j), where + ": MO/DM Hessians agree");
    }
}

chem::kumo::detail::DensityDerivatives from_mos(const chem::kumo::Molecule& mol,
                                                const Eigen::Vector3d& point,
                                                chem::kumo::SpinChannel spin) {
    chem::kumo::detail::DensityDerivatives out;
    const std::size_t occupied = mol.occupied_count(spin);
    for (std::size_t i = 1; i <= occupied; ++i) {
        const auto mo = kumo_test::Kernels::orbital_derivatives(mol, point, i, spin);
        out.value += mo.value * mo.value;
        out.gradient += 2.0 * mo.value * mo.gradient;
        out.hessian += 2.0 * (mo.gradient * mo.gradient.transpose() + mo.value * mo.hessian);
        out.tau += 0.5 * mo.gradient.squaredNorm();
    }
    out.laplacian = out.hessian.trace();
    return out;
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    return kumo_test::run([&] {
        for (int f = 1; f < argc; ++f) {
            const int failures_before = kumo_test::failures;
            const std::string name = argv[f];
            const auto mol = chem::kumo::Molecule::from_fchk(argv[f]);
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 13);
            for (int j = 0; j < points.cols(); ++j)
                points.col(j) = Eigen::Vector3d(0.31 + 0.21 * j, -0.47 + 0.17 * j, 0.39 - 0.11 * j);
            for (auto spin :
                 {chem::kumo::DensityComponent::Alpha, chem::kumo::DensityComponent::Beta}) {
                const auto block =
                    kumo_test::Kernels::density_derivatives_blocked(mol, points, spin, 5);
                const auto block1 =
                    kumo_test::Kernels::density_derivatives_blocked(mol, points, spin, 1);
                const std::string label = name + " spin " + std::to_string(static_cast<int>(spin));
                if (!kumo_test::check(block.size() == static_cast<std::size_t>(points.cols()),
                                      label + ": blocked output has one entry per point (got " +
                                          std::to_string(block.size()) + ")"))
                    continue;
                for (int j = 0; j < points.cols(); ++j) {
                    const auto r = points.col(j).eval();
                    const auto reference = kumo_test::Kernels::density_derivatives(mol, r, spin);
                    const std::string where = label + " point " + std::to_string(j);
                    compare(block[static_cast<std::size_t>(j)], reference, where + " block 5");
                    compare(block1[static_cast<std::size_t>(j)], reference, where + " block 1");
                    compare(from_mos(mol, r,
                                     spin == chem::kumo::DensityComponent::Alpha
                                         ? chem::kumo::SpinChannel::Alpha
                                         : chem::kumo::SpinChannel::Beta),
                            reference, where + " from MOs");
                }
            }
            for (auto spin :
                 {chem::kumo::DensityComponent::Total, chem::kumo::DensityComponent::Difference}) {
                const auto block =
                    kumo_test::Kernels::density_derivatives_blocked(mol, points, spin, 7);
                for (int j = 0; j < points.cols(); ++j)
                    compare(block[static_cast<std::size_t>(j)],
                            kumo_test::Kernels::density_derivatives(mol, points.col(j), spin),
                            name + " spin " + std::to_string(static_cast<int>(spin)) + " point " +
                                std::to_string(j) + " block 7");
            }
            kumo_test::check_throws<std::invalid_argument>(
                [&] {
                    (void)kumo_test::Kernels::density_derivatives_blocked(
                        mol, points, chem::kumo::DensityComponent::Total, 0);
                },
                name + ": density_derivatives_blocked rejects block size 0");
            if (kumo_test::failures == failures_before)
                std::cout << argv[f] << " MO/DM and blocked derivatives OK\n";
        }
    });
}
