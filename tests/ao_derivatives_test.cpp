#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/basis.hpp"
#include "chem/kumo/fchk.hpp"
#include "test_support.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    return kumo_test::run([&] {
        for (int file = 1; file < argc; ++file) {
            const std::string name = argv[file];
            const auto fchk = chem::kumo::FormattedCheckpoint(argv[file]);
            const chem::kumo::BasisEvaluator evaluator{chem::kumo::BasisSet(fchk)};
            double max_first = 0.0, max_second = 0.0, max_value = 0.0;
            constexpr double h = 1.0e-4;
            for (int t = 0; t < 12; ++t) {
                const Eigen::Vector3d r(0.31 + 0.37 * t, -0.24 + 0.19 * t, 0.13 - 0.21 * t);
                const auto exact = evaluator.evaluate_derivatives(r);
                max_value = std::max(max_value,
                                     (exact.values - evaluator.evaluate(r)).cwiseAbs().maxCoeff());
                for (int a = 0; a < 3; ++a) {
                    Eigen::Vector3d rp = r, rm = r;
                    rp[a] += h;
                    rm[a] -= h;
                    const auto plus = evaluator.evaluate_derivatives(rp);
                    const auto minus = evaluator.evaluate_derivatives(rm);
                    const auto fd1 = (plus.values - minus.values) / (2 * h);
                    max_first =
                        std::max(max_first, (fd1 - exact.gradient.col(a)).cwiseAbs().maxCoeff());
                    for (int b = a; b < 3; ++b) {
                        const int component = (a == b) ? a : (a == 0 ? (b == 1 ? 3 : 4) : 5);
                        const auto fd2 = (plus.gradient.col(b) - minus.gradient.col(b)) / (2 * h);
                        max_second = std::max(
                            max_second, (fd2 - exact.hessian.col(component)).cwiseAbs().maxCoeff());
                    }
                }
            }
            std::cout << argv[file] << " value=" << max_value << " grad=" << max_first
                      << " hessian=" << max_second << '\n';
            // The values from evaluate_derivatives() and evaluate() differ only
            // by rounding (and the 1e-14 AO cutoff): 1e-11 for O(1) AO values.
            kumo_test::check(!(max_value > 1.e-11),
                             name +
                                 ": AO values of evaluate_derivatives() match evaluate() "
                                 "within 1e-11 (max_abs=" +
                                 kumo_test::num(max_value) + ")");
            // Finite differences with h = 1e-4: central differences have
            // truncation error h^2/6 f''' (plus rounding ~eps |f| / h), and
            // tight, high-exponent primitives have large third derivatives,
            // hence the loose bounds (3e-3 gradient, 3e-2 Hessian). Wrong
            // derivative formulas give much larger errors.
            kumo_test::check(!(max_first > 3.e-3),
                             name +
                                 ": AO gradients match central differences within 3e-3 "
                                 "(max_abs=" +
                                 kumo_test::num(max_first) + ")");
            kumo_test::check(!(max_second > 3.e-2),
                             name +
                                 ": AO Hessians match central differences within 3e-2 "
                                 "(max_abs=" +
                                 kumo_test::num(max_second) + ")");
        }
    });
}
