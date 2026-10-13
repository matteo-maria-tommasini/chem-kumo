#include "chem/kumo/basis_evaluator.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    try {
        for (int a = 1; a < argc; ++a) {
            chem::kumo::FormattedCheckpoint fchk(argv[a]);
            chem::kumo::BasisEvaluator eval{chem::kumo::BasisSet(fchk)};
            Eigen::VectorXd fast(eval.basis_size()), reference(eval.basis_size());
            double worst = 0;
            for (int i = 0; i < 1000; ++i) {
                const double t = static_cast<double>(i);
                Eigen::Vector3d point(0.137*t - 10.0, std::sin(t*0.21)*4.0, std::cos(t*0.13)*3.0);
                eval.evaluate_into(point, fast);
                eval.evaluate_into_reference(point, reference);
                worst = std::max(worst, (fast-reference).cwiseAbs().maxCoeff());
            }
            std::cout << argv[a] << " basis=" << eval.basis_size() << " max_abs=" << worst << '\n';
            if (!std::isfinite(worst) || worst > 1e-12) return 1;
        }
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
