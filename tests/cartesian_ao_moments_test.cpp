#include "chem/kumo/basis_evaluator.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, 1, "FILE.fchk")) return 2;
    using kumo_test::check;
    using kumo_test::max_abs_diff;
    using kumo_test::num;
    return kumo_test::run([&] {
        chem::kumo::FormattedCheckpoint fchk(argv[1]);
        chem::kumo::BasisEvaluator ao{chem::kumo::BasisSet(fchk)};
        const auto m = ao.cartesian_moments();
        const double overlap_error = max_abs_diff(m.overlap, ao.overlap_matrix(), "overlap");
        check(overlap_error < 2e-12,
              "moment overlap matches the established overlap implementation within 2e-12 "
              "(max_abs=" +
                  num(overlap_error) + ")");
        const auto n = static_cast<Eigen::Index>(ao.basis_size());
        std::array<const Eigen::MatrixXd*, 10> matrices = {
            &m.overlap,   &m.first[0],  &m.first[1],  &m.first[2],  &m.second[0],
            &m.second[1], &m.second[2], &m.second[3], &m.second[4], &m.second[5]};
        for (std::size_t k = 0; k < matrices.size(); ++k) {
            const auto* matrix = matrices[k];
            const std::string label = "moment matrix " + std::to_string(k);
            if (!check(matrix->rows() == n && matrix->cols() == n,
                       label + " is " + std::to_string(n) + "x" + std::to_string(n) + " (got " +
                           std::to_string(matrix->rows()) + "x" + std::to_string(matrix->cols()) +
                           ")"))
                continue;
            check(matrix->allFinite(), label + " is finite");
            const double asymmetry = max_abs_diff(*matrix, matrix->transpose(), label);
            check(asymmetry < 2e-12,
                  label + " is symmetric within 2e-12 (max_abs=" + num(asymmetry) + ")");
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
                    check(std::abs(lhs - rhs) < 5e-10,
                          "AO " + std::to_string(i) + " axis " + std::to_string(d) +
                              ": first moment of a centered AO equals its center within 5e-10 "
                              "(expected " +
                              num(rhs, 17) + ", got " + num(lhs, 17) + ")");
                }
            }
            offset += count;
        }
        check(offset == n,
              "shells hold " + std::to_string(n) + " AOs (got " + std::to_string(offset) + ")");
        // For every AO, <(x-Rx)^2> must be nonnegative.
        offset = 0;
        for (const auto& shell : ao.basis().shells()) {
            const auto count = static_cast<Eigen::Index>(shell.function_count());
            for (Eigen::Index i = offset; i < offset + count; ++i) {
                for (int d = 0; d < 3; ++d) {
                    const double center = shell.origin[d];
                    const double variance =
                        m.second[static_cast<std::size_t>(d)](i, i) -
                        2 * center * m.first[static_cast<std::size_t>(d)](i, i) +
                        center * center * m.overlap(i, i);
                    check(variance > 0, "AO " + std::to_string(i) + " axis " + std::to_string(d) +
                                            ": spatial variance is positive (got " + num(variance) +
                                            ")");
                }
            }
            offset += count;
        }
        if (kumo_test::failures == 0)
            std::cout << "AO Cartesian moments validated for " << n << " functions\n";
    });
}
