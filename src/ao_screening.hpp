#pragma once
#include <Eigen/Core>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace chem::kumo::detail {
/// Check a screening threshold: finite and nonnegative (0 disables screening).
inline void validate_screening_threshold(double threshold) {
    if (!std::isfinite(threshold) || threshold < 0.0)
        throw std::invalid_argument("screening threshold must be finite and nonnegative");
}

/// Block screening of AO values. Gaussian AOs decay quickly, so in a block of
/// nearby grid points most of them are negligible everywhere (about 70% for
/// HBC def2-TZVP with 128-point blocks). Contracting only the significant AOs
/// replaces the nAO x nAO density-matrix product by a much smaller one.
///
/// Returns the rows of `ao` (one AO per row, one point per column) whose
/// absolute value exceeds `threshold` at least at one point of the block.
/// With threshold 0 every row is kept, so the contraction is exact.
///
/// Every dropped term of phi^T P phi contains a dropped AO factor of
/// magnitude <= threshold. With the default 1e-12 the observed change in the
/// density is at rounding level (about 1e-13 e/bohr^3 on HBC def2-TZVP, where
/// the density peaks at about 90 e/bohr^3).
inline std::vector<Eigen::Index> significant_aos(const Eigen::Ref<const Eigen::MatrixXd>& ao,
                                                 double threshold) {
    std::vector<Eigen::Index> rows;
    rows.reserve(static_cast<std::size_t>(ao.rows()));
    for (Eigen::Index row = 0; row < ao.rows(); ++row)
        if (threshold == 0.0 || ao.row(row).cwiseAbs().maxCoeff() > threshold) rows.push_back(row);
    return rows;
}

/// Copy of the given rows of `matrix`, in order. This and select_submatrix()
/// are written out rather than as Eigen indexed views, whose API differs
/// between Eigen 3.4 and Eigen 5 (for example `Eigen::all`).
inline Eigen::MatrixXd select_rows(const Eigen::Ref<const Eigen::MatrixXd>& matrix,
                                   const std::vector<Eigen::Index>& rows) {
    Eigen::MatrixXd result(static_cast<Eigen::Index>(rows.size()), matrix.cols());
    for (std::size_t i = 0; i < rows.size(); ++i)
        result.row(static_cast<Eigen::Index>(i)) = matrix.row(rows[i]);
    return result;
}

/// Copy of matrix(indices, indices): the rows and columns of `indices`, in order.
inline Eigen::MatrixXd select_submatrix(const Eigen::MatrixXd& matrix,
                                        const std::vector<Eigen::Index>& indices) {
    const auto n = static_cast<Eigen::Index>(indices.size());
    Eigen::MatrixXd result(n, n);
    for (Eigen::Index column = 0; column < n; ++column)
        for (Eigen::Index row = 0; row < n; ++row)
            result(row, column) = matrix(indices[static_cast<std::size_t>(row)],
                                         indices[static_cast<std::size_t>(column)]);
    return result;
}
} // namespace chem::kumo::detail
