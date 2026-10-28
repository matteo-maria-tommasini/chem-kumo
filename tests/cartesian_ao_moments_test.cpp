#include "chem/kumo/basis_evaluator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
double error(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
    return (a - b).cwiseAbs().maxCoeff();
}
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    try {
        chem::kumo::FormattedCheckpoint fchk(argv[1]);
        chem::kumo::BasisEvaluator ao{chem::kumo::BasisSet(fchk)};
        const auto m = ao.cartesian_moments();
        check(error(m.overlap, ao.overlap_matrix()) < 2e-12,
              "moment overlap differs from established overlap implementation");
        const auto n = static_cast<Eigen::Index>(ao.basis_size());
        std::array<const Eigen::MatrixXd*, 10> matrices = {
            &m.overlap, &m.first[0], &m.first[1], &m.first[2],
            &m.second[0], &m.second[1], &m.second[2],
            &m.second[3], &m.second[4], &m.second[5]};
        for (const auto* matrix : matrices) {
            check(matrix->rows() == n && matrix->cols() == n,
                  "unexpected moment matrix dimensions");
            check(matrix->allFinite(), "nonfinite moment matrix");
            check(error(*matrix, matrix->transpose()) < 2e-12,
                  "AO moment matrix is not symmetric");
        }

        // Every normalized Cartesian or pure AO centered at R has an even
        // probability density about R: its diagonal first moment is R.
        Eigen::Index offset = 0;
        for (const auto& shell : ao.basis().shells()) {
            const auto count = static_cast<Eigen::Index>(shell.function_count());
            for (Eigen::Index i = offset; i < offset + count; ++i) {
                for (int d = 0; d < 3; ++d) {
                    const double lhs = m.first[static_cast<std::size_t>(d)](i, i);
                    const double rhs = shell.origin[d] * m.overlap(i, i);
                    check(std::abs(lhs - rhs) < 5e-10, "incorrect first moment of centered AO");
                }
            }
            offset += count;
        }
        check(offset == n, "shell AO count mismatch");
        // For every AO, <(x-Rx)^2> must be nonnegative.
        offset = 0;
        for (const auto& shell : ao.basis().shells()) {
            const auto count = static_cast<Eigen::Index>(shell.function_count());
            for (Eigen::Index i = offset; i < offset + count; ++i) {
                for (int d = 0; d < 3; ++d) {
                    const double center = shell.origin[d];
                    const double variance = m.second[static_cast<std::size_t>(d)](i, i)
                        - 2 * center * m.first[static_cast<std::size_t>(d)](i, i)
                        + center * center * m.overlap(i, i);
                    check(variance > 0, "nonpositive AO spatial variance");
                }
            }
            offset += count;
        }
        std::cout << "AO Cartesian moments validated for " << n << " functions\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
