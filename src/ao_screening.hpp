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

/// Block screening of atomic-orbital (AO) values. Gaussian AOs decay quickly,
/// so in a block of nearby grid points most of them are negligible at every
/// point (about 70% for HBC def2-TZVP with 128-point blocks; HBC is the
/// 60-atom hexabenzocoronene test molecule in tests/data/hbc). Contracting
/// only the significant AOs replaces the nAO x nAO density-matrix product
/// (nAO = number of AOs) by a much smaller one.
///
/// Returns the rows of `ao` (one AO per row, one point per column) whose
/// absolute value exceeds `threshold` at least at one point of the block.
/// With threshold 0 every row is kept, so the contraction is exact. For
/// density derivatives, `ao` holds the value and derivative panels side by
/// side, so an AO is dropped only if its value and all its derivatives are
/// negligible; every dropped term of the derivative contractions then still
/// contains a factor of magnitude <= threshold.
///
/// In the density rho = phi^T P phi (phi = vector of AO values at a point,
/// P = AO density matrix), every dropped term contains a dropped AO factor of
/// magnitude <= threshold. With the default 1e-12 the observed change in the
/// density is at rounding level (about 1e-13 e/bohr^3 on HBC def2-TZVP, where
/// the density peaks at about 90 e/bohr^3).
inline std::vector<Eigen::Index> significant_aos(const Eigen::Ref<const Eigen::MatrixXd>& ao,
                                                 double threshold) {
    std::vector<Eigen::Index> rows;
    rows.reserve(static_cast<std::size_t>(ao.rows()));
    if (threshold == 0.0 || ao.cols() == 0) {
        for (Eigen::Index row = 0; row < ao.rows(); ++row)
            rows.push_back(row);
        return rows;
    }
    // Reduce column by column: `ao` is column-major, so this reads it in
    // memory order instead of striding along each row.
    const Eigen::VectorXd largest = ao.cwiseAbs().rowwise().maxCoeff();
    for (Eigen::Index row = 0; row < ao.rows(); ++row)
        if (largest[row] > threshold) rows.push_back(row);
    return rows;
}

/// Copy of the given rows of `matrix`, in order. This and select_submatrix()
/// are written out rather than as Eigen indexed views, whose API differs
/// between Eigen 3.4 and Eigen 5 (for example `Eigen::all`).
inline Eigen::MatrixXd select_rows(const Eigen::Ref<const Eigen::MatrixXd>& matrix,
                                   const std::vector<Eigen::Index>& rows) {
    const auto n = static_cast<Eigen::Index>(rows.size());
    Eigen::MatrixXd result(n, matrix.cols());
    // Column by column, in the memory order of both column-major matrices.
    for (Eigen::Index column = 0; column < matrix.cols(); ++column)
        for (Eigen::Index i = 0; i < n; ++i)
            result(i, column) = matrix(rows[static_cast<std::size_t>(i)], column);
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

/// The rows `kept` of `matrix`: a view of `matrix` itself when every row is
/// kept, otherwise a copy placed in `storage`.
inline Eigen::Ref<const Eigen::MatrixXd> kept_rows(const Eigen::Ref<const Eigen::MatrixXd>& matrix,
                                                   const std::vector<Eigen::Index>& kept,
                                                   Eigen::MatrixXd& storage) {
    if (static_cast<Eigen::Index>(kept.size()) == matrix.rows()) return matrix;
    storage = select_rows(matrix, kept);
    return storage;
}

/// matrix(kept, kept): `matrix` itself when every AO is kept, otherwise a
/// copy placed in `storage`.
inline const Eigen::MatrixXd& kept_submatrix(const Eigen::MatrixXd& matrix,
                                             const std::vector<Eigen::Index>& kept,
                                             Eigen::MatrixXd& storage) {
    if (static_cast<Eigen::Index>(kept.size()) == matrix.rows()) return matrix;
    storage = select_submatrix(matrix, kept);
    return storage;
}
} // namespace chem::kumo::detail
