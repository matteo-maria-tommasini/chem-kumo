#pragma once

namespace chem::kumo::detail {

/// Column of AODerivatives::hessian that holds the Cartesian second
/// derivative d2/(dr_a dr_b), for a, b in {0, 1, 2} = {x, y, z}.
///
/// The AO Hessian stores the six unique components as xx, yy, zz, xy, xz, yz
/// (see AODerivatives in include/chem/kumo/basis_evaluator.hpp). Note that the Evaluator output
/// uses a different, row-major order (xx, xy, xz, yy, yz, zz); do not mix the two.
inline constexpr int ao_hessian_column[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};

} // namespace chem::kumo::detail
