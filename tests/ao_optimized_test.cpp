#include "chem/kumo/basis_evaluator.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    return kumo_test::run([&] {
        for (int a = 1; a < argc; ++a) {
            chem::kumo::FormattedCheckpoint fchk(argv[a]);
            chem::kumo::BasisEvaluator eval{chem::kumo::BasisSet(fchk)};
            Eigen::VectorXd fast(eval.basis_size()), reference(eval.basis_size());
            double worst = 0;
            for (int i = 0; i < 1000; ++i) {
                const double t = static_cast<double>(i);
                Eigen::Vector3d point(0.137 * t - 10.0, std::sin(t * 0.21) * 4.0,
                                      std::cos(t * 0.13) * 3.0);
                eval.evaluate_into(point, fast);
                eval.evaluate_into_reference(point, reference);
                worst = std::max(worst, kumo_test::max_abs_diff(fast, reference));
            }
            std::cout << argv[a] << " basis=" << eval.basis_size() << " max_abs=" << worst << '\n';
            // The optimized kernel differs from the original pointwise
            // evaluator only by rounding and by zeroing AOs below
            // BasisEvaluator::ao_cutoff (1e-14). AO values are O(1), so 1e-12
            // allows rounding differences of a few thousand ulps.
            kumo_test::check(std::isfinite(worst) && worst <= 1e-12,
                             std::string(argv[a]) +
                                 ": optimized AO values match the reference kernel within 1e-12 "
                                 "(max_abs=" +
                                 kumo_test::num(worst) + ")");
        }
    });
}
