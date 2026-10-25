#pragma once

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chem::kumo {

/// S-orthogonal projector onto the column space of an AO coefficient matrix.
/// Stores an S-orthonormal basis, not an M-by-M projector.
/// The AO overlap metric must be symmetric positive definite.
class SubspaceProjector {
public:
    SubspaceProjector(const Eigen::Ref<const Eigen::MatrixXd>& metric,
                      const Eigen::Ref<const Eigen::MatrixXd>& columns,
                      double relative_rank_tolerance = 1e-12)
        : metric_(metric) {
        if (metric_.rows() == 0 || metric_.rows() != metric_.cols() ||
            columns.rows() != metric_.rows() || columns.cols() == 0 ||
            columns.cols() > columns.rows() ||
            !metric_.allFinite() || !columns.allFinite() ||
            !std::isfinite(relative_rank_tolerance) ||
            relative_rank_tolerance <= 0.0 || relative_rank_tolerance >= 1.0)
            throw std::invalid_argument("SubspaceProjector: invalid dimensions, data, or tolerance");

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
            throw std::invalid_argument("SubspaceProjector: rank-deficient or ill-conditioned subspace");

        // Symmetric inverse square root avoids an explicit matrix inverse.
        const Eigen::MatrixXd inverse_sqrt =
            eig.eigenvectors() * values.array().rsqrt().matrix().asDiagonal()
            * eig.eigenvectors().transpose();
        basis_.noalias() = columns * inverse_sqrt;
    }

    [[nodiscard]] Eigen::Index ambient_dimension() const noexcept { return metric_.rows(); }
    [[nodiscard]] Eigen::Index dimension() const noexcept { return basis_.cols(); }
    [[nodiscard]] const Eigen::MatrixXd& orthonormal_basis() const noexcept { return basis_; }

    /// Apply to one or several columns of AO expansion coefficients.
    [[nodiscard]] Eigen::MatrixXd apply(const Eigen::Ref<const Eigen::MatrixXd>& vectors) const {
        if (vectors.rows() != metric_.rows())
            throw std::invalid_argument("SubspaceProjector: vector dimension mismatch");
        if (!vectors.allFinite())
            throw std::invalid_argument("SubspaceProjector: non-finite projection input");
        return basis_ * (basis_.transpose() * (metric_ * vectors));
    }

    /// Materialize the coefficient-space projector, for tests and diagnostics.
    [[nodiscard]] Eigen::MatrixXd matrix() const {
        return basis_ * (basis_.transpose() * metric_);
    }

private:
    Eigen::MatrixXd metric_;
    Eigen::MatrixXd basis_;
};

} // namespace chem::kumo
