#pragma once

#include "chem/kumo/basis_evaluator.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace chem::kumo {

/// Accumulates the atomic-orbital (AO) matrix of a scalar field f by a
/// quadrature supplied by the caller,
///   W_mu,nu = sum_q w_q f(r_q) phi_mu(r_q) phi_nu(r_q),
/// an approximation of the integral of f phi_mu phi_nu over space. The
/// points r_q are in bohr, the quadrature weights w_q include the volume
/// element (bohr^3), and f(r_q) are the field values; neither weights nor
/// values are normalized here. The caller chooses the points and the field.
class WeightedOverlapAccumulator {
public:
    explicit WeightedOverlapAccumulator(const BasisEvaluator& evaluator)
        : evaluator_(evaluator),
          matrix_(Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(evaluator.basis_size()),
                                        static_cast<Eigen::Index>(evaluator.basis_size()))) {}

    /// Add the contribution of one block of points: `points` is 3 x n, and
    /// `quadrature_weights` and `field_values` have n entries. An empty block
    /// adds nothing. Throws std::invalid_argument for mismatched sizes or a
    /// nonfinite entry.
    void add(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
             const Eigen::Ref<const Eigen::VectorXd>& quadrature_weights,
             const Eigen::Ref<const Eigen::VectorXd>& field_values) {
        const Eigen::Index n = points.cols();
        if (quadrature_weights.size() != n || field_values.size() != n || !points.allFinite() ||
            !quadrature_weights.allFinite() || !field_values.allFinite())
            throw std::invalid_argument("WeightedOverlapAccumulator: invalid quadrature data");
        if (n == 0) return;
        Eigen::MatrixXd ao(static_cast<Eigen::Index>(evaluator_.basis_size()), n);
        evaluator_.evaluate_batch(points, ao);
        const Eigen::VectorXd effective =
            (quadrature_weights.array() * field_values.array()).matrix();
        // Scale the columns of a copy of the AO block by w_q f(r_q) instead
        // of multiplying by an n x n diagonal matrix.
        Eigen::MatrixXd scaled = ao;
        scaled.array().rowwise() *= effective.transpose().array();
        matrix_.noalias() += scaled * ao.transpose();
    }

    /// Accumulated nAO x nAO matrix.
    [[nodiscard]] const Eigen::MatrixXd& matrix() const noexcept { return matrix_; }

    /// Clear the accumulated matrix.

    void reset() noexcept { matrix_.setZero(); }

private:
    const BasisEvaluator& evaluator_; ///< must outlive the accumulator
    Eigen::MatrixXd matrix_;
};

/// The matrix of WeightedOverlapAccumulator for all points at once, added
/// in blocks of `block_size` points to bound the memory used. Throws
/// std::invalid_argument for a zero block size, mismatched sizes or a
/// nonfinite entry.
[[nodiscard]] inline Eigen::MatrixXd
weighted_overlap_matrix(const BasisEvaluator& evaluator,
                        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                        const Eigen::Ref<const Eigen::VectorXd>& quadrature_weights,
                        const Eigen::Ref<const Eigen::VectorXd>& field_values,
                        std::size_t block_size = 128) {
    if (block_size == 0) throw std::invalid_argument("weighted_overlap_matrix: zero block size");
    if (quadrature_weights.size() != points.cols() || field_values.size() != points.cols())
        throw std::invalid_argument("weighted_overlap_matrix: inconsistent quadrature sizes");
    WeightedOverlapAccumulator accumulator(evaluator);
    const Eigen::Index chunk = static_cast<Eigen::Index>(
        std::min(block_size, static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max())));
    for (Eigen::Index first = 0; first < points.cols();) {
        const Eigen::Index count = std::min(chunk, points.cols() - first);
        accumulator.add(points.middleCols(first, count), quadrature_weights.segment(first, count),
                        field_values.segment(first, count));
        first += count;
    }
    return accumulator.matrix();
}

} // namespace chem::kumo
