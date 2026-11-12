// BasisEvaluator::evaluate_derivatives_batch: every panel of a batch agrees
// with the single-point evaluate_derivatives() at that point (orders 1 and
// 2), the value panel agrees with evaluate(), and invalid orders or output
// shapes are rejected.
#include "chem/kumo/basis.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/fchk.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
// Largest |a - b|, relative to the largest |a| in the column (or 1 if smaller).
double relative_difference(const Eigen::Ref<const Eigen::VectorXd>& a,
                           const Eigen::Ref<const Eigen::VectorXd>& b) {
    return (a - b).cwiseAbs().maxCoeff() / std::max(1.0, a.cwiseAbs().maxCoeff());
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " FILE.fchk...\n";
        return 2;
    }
    using chem::kumo::BasisEvaluator;
    for (int file = 1; file < argc; ++file) {
        const std::string name = argv[file];
        const BasisEvaluator evaluator{chem::kumo::BasisSet(chem::kumo::FormattedCheckpoint(name))};
        const auto n = static_cast<Eigen::Index>(evaluator.basis_size());
        // 37 points spread over a few bohr, including one at the origin.
        constexpr Eigen::Index count = 37;
        Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, count);
        for (Eigen::Index j = 0; j < count; ++j)
            points.col(j) << 0.29 * j - 4.0, 2.5 * std::sin(0.7 * j), 1.9 * std::cos(1.3 * j);
        points.col(5).setZero();

        for (int order : {1, 2}) {
            const int panels = BasisEvaluator::derivative_panel_count(order);
            kumo_test::check(panels == (order == 2 ? 10 : 4), "derivative_panel_count(" +
                                                                  std::to_string(order) + ") is " +
                                                                  (order == 2 ? "10" : "4"));
            Eigen::MatrixXd batch(n, panels * count);
            evaluator.evaluate_derivatives_batch(points, order, batch);
            double error = 0.0, value_error = 0.0;
            for (Eigen::Index j = 0; j < count; ++j) {
                const Eigen::Vector3d point = points.col(j);
                const auto single = evaluator.evaluate_derivatives(point, order);
                const auto panel = [&](int c) { return batch.col(c * count + j); };
                error = std::max(error, relative_difference(single.values, panel(0)));
                for (int a = 0; a < 3; ++a)
                    error =
                        std::max(error, relative_difference(single.gradient.col(a), panel(1 + a)));
                for (int h = 0; order == 2 && h < 6; ++h)
                    error =
                        std::max(error, relative_difference(single.hessian.col(h), panel(4 + h)));
                value_error =
                    std::max(value_error, relative_difference(evaluator.evaluate(point), panel(0)));
            }
            const std::string label = name + " order " + std::to_string(order);
            kumo_test::check(error <= 1e-14, label +
                                                 ": batch panels match evaluate_derivatives() "
                                                 "within 1e-14 relative (max " +
                                                 kumo_test::num(error) + ")");
            kumo_test::check(value_error <= 1e-14, label +
                                                       ": value panel matches evaluate() within "
                                                       "1e-14 relative (max " +
                                                       kumo_test::num(value_error) + ")");
        }

        Eigen::MatrixXd wrong(n, 4 * count - 1);
        kumo_test::check_throws<std::runtime_error>(
            [&] { evaluator.evaluate_derivatives_batch(points, 1, wrong); },
            name + ": wrong output width is rejected");
        Eigen::MatrixXd short_rows(n - 1, 4 * count);
        kumo_test::check_throws<std::runtime_error>(
            [&] { evaluator.evaluate_derivatives_batch(points, 1, short_rows); },
            name + ": wrong output height is rejected");
        for (int order : {0, 3}) {
            Eigen::MatrixXd output(n, 10 * count);
            kumo_test::check_throws<std::invalid_argument>(
                [&] { evaluator.evaluate_derivatives_batch(points, order, output); },
                name + ": derivative order " + std::to_string(order) + " is rejected");
        }
        Eigen::Matrix<double, 3, Eigen::Dynamic> none(3, 0);
        Eigen::MatrixXd empty(n, 0);
        evaluator.evaluate_derivatives_batch(none, 2, empty);
    }
    return kumo_test::exit_code();
}
