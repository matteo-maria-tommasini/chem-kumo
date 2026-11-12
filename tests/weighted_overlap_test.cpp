#include "chem/kumo/weighted_overlap.hpp"
#include "chem/kumo/fchk.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
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
            // Pass if error < bound; the message shows both.
            const auto below = [&](double error, double bound, const std::string& what) {
                kumo_test::check(error < bound, name + what + " (max_abs=" + num(error) +
                                                    ", tolerance " + num(bound) + ")");
            };
            FormattedCheckpoint fchk(argv[arg]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            constexpr Eigen::Index n = 17;
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, n);
            Eigen::VectorXd weights(n), field(n), ones = Eigen::VectorXd::Ones(n);
            for (Eigen::Index j = 0; j < n; ++j) {
                points.col(j) = Eigen::Vector3d(0.13 * j - 0.8, 0.17 * j - 0.5, 0.09 * j + 0.2);
                weights[j] = 0.01 + 0.002 * j;
                field[j] = 0.3 + 0.07 * j;
            }
            Eigen::MatrixXd ao(evaluator.basis_size(), n);
            evaluator.evaluate_batch(points, ao);
            const Eigen::VectorXd effective = (weights.array() * field.array()).matrix();
            const Eigen::MatrixXd reference = ao * effective.asDiagonal() * ao.transpose();
            const Eigen::MatrixXd blocked =
                weighted_overlap_matrix(evaluator, points, weights, field, 5);
            const Eigen::MatrixXd singleton =
                weighted_overlap_matrix(evaluator, points, weights, field, 1);
            const double tolerance = 3e-12 * std::max(1.0, reference.cwiseAbs().maxCoeff());
            below(max_abs_diff(reference, blocked), tolerance,
                  "blocked matrix (block 5) matches the dense reference");
            below(max_abs_diff(reference, singleton), tolerance,
                  "single-point blocks match the dense reference");
            below(max_abs_diff(blocked, blocked.transpose()), tolerance,
                  "blocked matrix is symmetric");
            const Eigen::VectorXd double_field = 2.0 * field;
            const Eigen::MatrixXd twice =
                weighted_overlap_matrix(evaluator, points, weights, double_field, 7);
            below(max_abs_diff(twice, 2.0 * reference), 2 * tolerance,
                  "doubling the field doubles the matrix");
            const Eigen::VectorXd zero_field = Eigen::VectorXd::Zero(n);
            const Eigen::MatrixXd zero =
                weighted_overlap_matrix(evaluator, points, weights, zero_field, 7);
            kumo_test::check(zero.cwiseAbs().maxCoeff() == 0.0,
                             name + "a zero field gives exactly the zero matrix (max_abs=" +
                                 num(zero.cwiseAbs().maxCoeff()) + ")");
            // Positive semidefiniteness without a costly full eigendecomposition.
            for (int k = 0; k < 4; ++k) {
                Eigen::VectorXd v = Eigen::VectorXd::Zero(evaluator.basis_size());
                for (Eigen::Index i = 0; i < v.size(); ++i)
                    v[i] = std::sin(0.17 * (i + 1) * (k + 1));
                const double quadratic = v.dot(blocked * v);
                const double bound = -tolerance * v.squaredNorm();
                kumo_test::check(quadratic >= bound, name + "v^T M v >= " + num(bound) +
                                                         " for test vector " + std::to_string(k) +
                                                         " (positive semidefinite; got " +
                                                         num(quadratic) + ")");
            }
            // Signed weights are permitted and must not be silently clipped.
            const Eigen::VectorXd negative_field = -field;
            const Eigen::MatrixXd negative =
                weighted_overlap_matrix(evaluator, points, weights, negative_field, 6);
            below(max_abs_diff(negative, -reference), tolerance,
                  "a negated field negates the matrix");
            WeightedOverlapAccumulator accumulator(evaluator);
            accumulator.add(points, weights, field);
            below(max_abs_diff(accumulator.matrix(), reference), tolerance,
                  "WeightedOverlapAccumulator matches the dense reference");
            accumulator.reset();
            kumo_test::check(accumulator.matrix().cwiseAbs().maxCoeff() == 0.0,
                             name + "reset() clears the accumulator exactly");
            kumo_test::check_throws<std::invalid_argument>(
                [&] { (void)weighted_overlap_matrix(evaluator, points, weights, field, 0); },
                name + "weighted_overlap_matrix rejects block size 0");
            std::cout << argv[arg] << " nAO=" << evaluator.basis_size() << " points=" << n
                      << " max_abs=" << max_abs_diff(reference, blocked) << '\n';
        }
    });
}
