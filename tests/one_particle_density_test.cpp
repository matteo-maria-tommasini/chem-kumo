#include "chem/kumo/density.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    using kumo_test::max_abs_diff;
    using kumo_test::num;
    return kumo_test::run([&] {
        for (int arg = 1; arg < argc; ++arg) {
            const std::string name = std::string(argv[arg]) + ": ";
            const auto m = Molecule::from_fchk(argv[arg]);
            FormattedCheckpoint f(argv[arg]);
            const Eigen::MatrixXd S = BasisEvaluator{BasisSet(f)}.overlap_matrix();
            const auto pt = m.density_matrix(Spin::Total);
            const auto ps = m.density_matrix(Spin::Difference);
            const auto pa = m.density_matrix(Spin::Alpha);
            const auto pb = m.density_matrix(Spin::Beta);
            const double total_split = max_abs_diff(pa + pb, pt);
            const double spin_split = max_abs_diff(pa - pb, ps);
            kumo_test::check(!(total_split > 2e-14),
                             name + "P_total = P_alpha + P_beta within 2e-14 (max_abs=" +
                                 num(total_split) + ")");
            kumo_test::check(
                !(spin_split > 2e-14),
                name + "P_spin = P_alpha - P_beta within 2e-14 (max_abs=" + num(spin_split) + ")");
            const auto verify = [&](const Eigen::MatrixXd& p, double electrons, const char* label) {
                const std::string where = name + label + " 1-RDM: ";
                const double asymmetry = max_abs_diff(p, p.transpose());
                kumo_test::check(
                    p.allFinite() && !(asymmetry > 1e-12),
                    where + "finite and symmetric within 1e-12 (max_abs=" + num(asymmetry) + ")");
                const double trace = (p.cwiseProduct(S)).sum();
                const double idempotency = max_abs_diff(p * S * p, p);
                // The SCF density is read from FCHK; it is not overwritten
                // by an occupied-orbital reconstruction. Tolerances account
                // for printed FCHK precision and AO metric conditioning.
                kumo_test::check(!(std::abs(trace - electrons) > 2e-4),
                                 where + "tr(P S) equals the electron count " + num(electrons) +
                                     " within 2e-4 (got " + num(trace, 12) + ")");
                kumo_test::check(
                    !(m.electronic_structure_info().ao_mo_eligible() && idempotency > 2e-4),
                    where + "P S P = P within 2e-4 (max_abs=" + num(idempotency) + ")");
                std::cout << argv[arg] << " spin=" << label
                          << " electrons=" << std::setprecision(12) << trace
                          << " idempotency_max=" << idempotency << '\n';
            };
            verify(pa, m.occupied_count(Spin::Alpha), "alpha");
            verify(pb, m.occupied_count(Spin::Beta), "beta");
            const Eigen::Vector3d point(0.29, -0.41, 0.37);
            const Eigen::VectorXd ao = BasisEvaluator{BasisSet(f)}.evaluate(point);
            for (Spin spin : {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference}) {
                const Eigen::MatrixXd p = m.density_matrix(spin);
                const double direct = ao.dot(p * ao);
                const double reference = m.density_derivatives(point, spin).value;
                const double tolerance =
                    1e-9 * std::max({1.0, std::abs(direct), std::abs(reference)});
                kumo_test::check(!(std::abs(direct - reference) > tolerance),
                                 name + "spin " + std::to_string(static_cast<int>(spin)) +
                                     ": phi^T P phi matches density_derivatives() (expected " +
                                     num(reference, 17) + ", got " + num(direct, 17) +
                                     ", tolerance " + num(tolerance) + ")");
            }
        }
    });
}
