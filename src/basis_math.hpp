#pragma once

#include "chem/kumo/basis.hpp"
#include "chem/kumo/basis_evaluator.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chem::kumo::detail {

/// n!! = n (n - 2) (n - 4) ... for odd n; 1 for n <= 0 (so (-1)!! = 1).
double double_factorial_odd(int n) noexcept;
/// Normalization constant of the Cartesian Gaussian primitive
/// x^i y^j z^k exp(-alpha r^2) with (i, j, k) = `power`.
double cartesian_primitive_normalization(double alpha, const CartesianPower& power);

/// AOs of one angular momentum on one centre, with normalized contraction
/// weights. SP shells give an s block followed by a p block, so the blocks
/// cover the AOs in FCHK order.
struct NormalizedBlock {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    int angular_momentum = 0;
    std::vector<double> exponents;
    std::vector<CartesianPower> powers; ///< Cartesian components, Gaussian order
    /// Contraction coefficient times primitive normalization; one row per
    /// Cartesian component, one column per primitive.
    Eigen::MatrixXd weights;
    /// Cartesian-to-pure (spherical) transformation, one row per pure AO and
    /// one column per Cartesian component; empty for Cartesian blocks and for
    /// s and p blocks, where both representations coincide.
    Eigen::MatrixXd pure_transform;
    std::size_t offset = 0; ///< index of the block's first AO in the basis
    std::size_t count = 0;  ///< number of AOs in the block
};

/// The blocks of a basis; throws std::logic_error if they do not cover
/// exactly basis.function_count() AOs.
std::vector<NormalizedBlock> normalized_blocks(const BasisSet& basis);

} // namespace chem::kumo::detail
