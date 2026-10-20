#pragma once

#include "chem/kumo/basis_evaluator.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <limits>

namespace chem::kumo {

/// Accumulate an AO weighted-overlap matrix from a supplied spatial quadrature.
/// Points are 3 x N (bohr); quadrature_weights include volume elements (bohr^3);
/// field_values are sampled scalar multipliers. Neither array is normalized here.
/// The result is a quadrature approximation, not an analytic integral.
/// This class deliberately does not choose a spatial grid or interpret the field.
class WeightedOverlapAccumulator {
public:
    explicit WeightedOverlapAccumulator(const BasisEvaluator& evaluator)
        : evaluator_(evaluator), matrix_(Eigen::MatrixXd::Zero(
              static_cast<Eigen::Index>(evaluator.basis_size()),
              static_cast<Eigen::Index>(evaluator.basis_size()))) {}

    /// Add one point block. Empty blocks are allowed.
    void add(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
             const Eigen::Ref<const Eigen::VectorXd>& quadrature_weights,
             const Eigen::Ref<const Eigen::VectorXd>& field_values) {
        const Eigen::Index n = points.cols();
        if (quadrature_weights.size() != n || field_values.size() != n ||
            !points.allFinite() || !quadrature_weights.allFinite() ||
            !field_values.allFinite())
            throw std::invalid_argument("WeightedOverlapAccumulator: invalid quadrature data");
        if (n == 0) return;
        Eigen::MatrixXd ao(static_cast<Eigen::Index>(evaluator_.basis_size()), n);
        evaluator_.evaluate_batch(points, ao);
        const Eigen::VectorXd effective =
            (quadrature_weights.array() * field_values.array()).matrix();
        // Scale only the temporary AO block. Do not build an N x N diagonal matrix.
        Eigen::MatrixXd scaled = ao;
        scaled.array().rowwise() *= effective.transpose().array();
        matrix_.noalias() += scaled * ao.transpose();
    }

    [[nodiscard]] const Eigen::MatrixXd& matrix() const noexcept { return matrix_; }
    void reset() noexcept { matrix_.setZero(); }

private:
    const BasisEvaluator& evaluator_; // evaluator must outlive the accumulator
    Eigen::MatrixXd matrix_;
};

/// Convenience blocked integration for an externally supplied quadrature.
[[nodiscard]] inline Eigen::MatrixXd weighted_overlap_matrix(
    const BasisEvaluator& evaluator,
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
        accumulator.add(points.middleCols(first, count),
                        quadrature_weights.segment(first, count),
                        field_values.segment(first, count));
        first += count;
    }
    return accumulator.matrix();
}

} // namespace chem::kumo
