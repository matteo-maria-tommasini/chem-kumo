#pragma once

#include "chem/kumo/basis.hpp"
#include "chem/kumo/basis_evaluator.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chem::kumo::detail {

double double_factorial_odd(int n) noexcept;
double cartesian_primitive_normalization(double alpha, const CartesianPower& power);

/// AOs of one angular momentum on one centre, with normalized contraction
/// weights. SP shells give an s block followed by a p block, so the blocks
/// cover the AOs in FCHK order.
struct NormalizedBlock {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    int angular_momentum = 0;
    std::vector<double> exponents;
    std::vector<CartesianPower> powers; ///< Cartesian components, Gaussian order
    Eigen::MatrixXd weights;            ///< Cartesian component x primitive
    Eigen::MatrixXd pure_transform;     ///< pure x Cartesian, empty if Cartesian
    std::size_t offset = 0;             ///< first AO
    std::size_t count = 0;              ///< number of AOs
};

/// The blocks of a basis; throws std::logic_error if they do not cover
/// exactly basis.function_count() AOs.
std::vector<NormalizedBlock> normalized_blocks(const BasisSet& basis);

} // namespace chem::kumo::detail
