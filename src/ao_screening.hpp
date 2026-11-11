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
} // namespace chem::kumo::detail
