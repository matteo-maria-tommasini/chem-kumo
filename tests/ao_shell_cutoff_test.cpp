// Shell cutoff of the AO evaluators: wherever evaluate() returns an exact zero
// for an AO that the reference evaluator (no cutoff) does not, the true value
// is at most BasisEvaluator::ao_cutoff, finite differences of the reference
// show that the gradient and the diagonal second derivatives are as small,
// and evaluate_derivatives_batch() returns exact zeros in every panel. The
// test also checks that the cutoff is exercised.
#include "chem/kumo/basis.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/fchk.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        using chem::kumo::BasisEvaluator;
        const double cutoff = BasisEvaluator::ao_cutoff;
        for (int file = 1; file < argc; ++file) {
            const std::string name = argv[file];
            const BasisEvaluator evaluator{
                chem::kumo::BasisSet(chem::kumo::FormattedCheckpoint(name))};
            const auto n = static_cast<Eigen::Index>(evaluator.basis_size());
            // Points along three directions out to 40 bohr, so that every shell
            // is cut somewhere and many are cut close to their radius.
            constexpr Eigen::Index count = 240;
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, count);
            const Eigen::Vector3d directions[3] = {Eigen::Vector3d(1.0, 0.3, -0.2).normalized(),
                                                   Eigen::Vector3d(-0.4, 1.0, 0.5).normalized(),
                                                   Eigen::Vector3d(0.2, -0.6, 1.0).normalized()};
            for (Eigen::Index j = 0; j < count; ++j)
                points.col(j) = directions[j % 3] * (0.5 * static_cast<double>(j / 3));

            Eigen::MatrixXd panels(n, 10 * count);
            evaluator.evaluate_derivatives_batch(points, 2, panels);
            constexpr double h = 1e-3;
            long cut = 0;
            double value_bound = 0.0, gradient_bound = 0.0, second_bound = 0.0;
            bool derivative_zeros = true;
            Eigen::VectorXd fast(n), reference(n), plus(n), minus(n);
            for (Eigen::Index j = 0; j < count; ++j) {
                const Eigen::Vector3d point = points.col(j);
                evaluator.evaluate_into(point, fast);
                evaluator.evaluate_into_reference(point, reference);
                for (Eigen::Index mu = 0; mu < n; ++mu) {
                    if (fast[mu] != 0.0 || reference[mu] == 0.0) continue;
                    ++cut;
                    value_bound = std::max(value_bound, std::abs(reference[mu]));
                    for (int c = 0; c < 10; ++c)
                        derivative_zeros = derivative_zeros && panels(mu, c * count + j) == 0.0;
                    for (int axis = 0; axis < 3; ++axis) {
                        Eigen::Vector3d step = Eigen::Vector3d::Zero();
                        step[axis] = h;
                        evaluator.evaluate_into_reference(point + step, plus);
                        evaluator.evaluate_into_reference(point - step, minus);
                        gradient_bound =
                            std::max(gradient_bound, std::abs(plus[mu] - minus[mu]) / (2.0 * h));
                        second_bound = std::max(
                            second_bound,
                            std::abs(plus[mu] - 2.0 * reference[mu] + minus[mu]) / (h * h));
                    }
                }
            }
            kumo_test::check(cut > 0, name + ": the shell cutoff zeroes some nonzero AOs");
            kumo_test::check(value_bound <= cutoff,
                             name + ": cut AO values are at most ao_cutoff (max " +
                                 kumo_test::num(value_bound) + ")");
            // Finite differences add a small truncation error, hence the factor 2.
            kumo_test::check(gradient_bound <= 2.0 * cutoff,
                             name + ": cut AO gradients are at most ao_cutoff (max " +
                                 kumo_test::num(gradient_bound) + ")");
            kumo_test::check(second_bound <= 2.0 * cutoff,
                             name + ": cut AO second derivatives are at most ao_cutoff (max " +
                                 kumo_test::num(second_bound) + ")");
            kumo_test::check(derivative_zeros,
                             name + ": evaluate_derivatives_batch() returns zeros for cut AOs");
            std::cout << name << ": " << cut << " cut AO values; max value " << value_bound
                      << ", gradient " << gradient_bound << ", second derivative " << second_bound
                      << "\n";
        }
    });
}
