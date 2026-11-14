#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chem::kumo {

/// Projector onto the space spanned by the columns of an atomic-orbital (AO)
/// coefficient matrix C, such as a set of molecular orbitals, orthogonal
/// with respect to the AO overlap matrix S (the metric):
///   P v = B B^T S v,  with B = C (C^T S C)^(-1/2),
/// so that B^T S B = 1 (B is an S-orthonormal basis of the space). Only B is
/// stored, not the M x M matrix of P (M = number of AOs).
class SubspaceProjector {
public:
    /// Throws std::invalid_argument for empty or mismatched dimensions, more
    /// columns than AOs, nonfinite data, a tolerance outside (0, 1), a metric
    /// that is not symmetric (to 1e-12 relative) or not positive definite,
    /// or columns that are not linearly independent: an eigenvalue of
    /// C^T S C at most relative_rank_tolerance times the largest one. Throws
    /// std::runtime_error if the eigenvalue solver fails.
    SubspaceProjector(const Eigen::Ref<const Eigen::MatrixXd>& metric,
                      const Eigen::Ref<const Eigen::MatrixXd>& columns,
                      double relative_rank_tolerance = 1e-12)
        : metric_(metric) {
        if (metric_.rows() == 0 || metric_.rows() != metric_.cols() ||
            columns.rows() != metric_.rows() || columns.cols() == 0 ||
            columns.cols() > columns.rows() || !metric_.allFinite() || !columns.allFinite() ||
            !std::isfinite(relative_rank_tolerance) || relative_rank_tolerance <= 0.0 ||
            relative_rank_tolerance >= 1.0)
            throw std::invalid_argument(
                "SubspaceProjector: invalid dimensions, data, or tolerance");

        const double symmetry_error = (metric_ - metric_.transpose()).cwiseAbs().maxCoeff();
        if (symmetry_error > 1e-12 * std::max(1.0, metric_.cwiseAbs().maxCoeff()))
            throw std::invalid_argument("SubspaceProjector: metric is not symmetric");
        metric_ = 0.5 * (metric_ + Eigen::MatrixXd(metric_.transpose()));
        Eigen::LLT<Eigen::MatrixXd> metric_llt(metric_);
        if (metric_llt.info() != Eigen::Success)
            throw std::invalid_argument("SubspaceProjector: metric is not positive definite");

        const Eigen::MatrixXd gram = columns.transpose() * metric_ * columns;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(gram);
        if (eig.info() != Eigen::Success)
            throw std::runtime_error("SubspaceProjector: Gram eigensolver failed");
        const auto& values = eig.eigenvalues();
        const double largest = values.maxCoeff();
        if (!(largest > 0.0) || values.minCoeff() <= relative_rank_tolerance * largest)
            throw std::invalid_argument(
                "SubspaceProjector: rank-deficient or ill-conditioned subspace");

        // (C^T S C)^(-1/2) from the eigendecomposition of C^T S C, without
        // an explicit matrix inverse.
        const Eigen::MatrixXd inverse_sqrt = eig.eigenvectors() *
                                             values.array().rsqrt().matrix().asDiagonal() *
                                             eig.eigenvectors().transpose();
        basis_.noalias() = columns * inverse_sqrt;
    }

    /// Number of AOs.
    [[nodiscard]] Eigen::Index ambient_dimension() const noexcept { return metric_.rows(); }

    /// Dimension of the subspace: the number of input columns, which must be
    /// linearly independent within the rank tolerance (otherwise the
    /// constructor throws).
    [[nodiscard]] Eigen::Index dimension() const noexcept { return basis_.cols(); }

    /// S-orthonormal basis of the subspace, one AO coefficient vector per column.
    [[nodiscard]] const Eigen::MatrixXd& orthonormal_basis() const noexcept { return basis_; }

    /// Projection P v of each column v of AO coefficients in `vectors`.
    /// Throws std::invalid_argument if the number of rows differs from the
    /// number of AOs or an entry is not finite.
    [[nodiscard]] Eigen::MatrixXd apply(const Eigen::Ref<const Eigen::MatrixXd>& vectors) const {
        if (vectors.rows() != metric_.rows())
            throw std::invalid_argument("SubspaceProjector: vector dimension mismatch");
        if (!vectors.allFinite())
            throw std::invalid_argument("SubspaceProjector: non-finite projection input");
        return basis_ * (basis_.transpose() * (metric_ * vectors));
    }

    /// The M x M matrix of the projector, B B^T S, for tests and diagnostics.
    [[nodiscard]] Eigen::MatrixXd matrix() const { return basis_ * (basis_.transpose() * metric_); }

private:
    Eigen::MatrixXd metric_;
    Eigen::MatrixXd basis_;
};

} // namespace chem::kumo
