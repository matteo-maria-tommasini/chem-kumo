#pragma once

namespace chem::kumo::detail {

/// Column of AODerivatives::hessian (and offset of the second-derivative
/// panels of BasisEvaluator::evaluate_derivatives_batch(), after the value and
/// gradient panels) that holds the Cartesian second derivative d2/(dr_a dr_b),
/// for a, b in {0, 1, 2} = {x, y, z}.
///
/// Every six-component Hessian of the public API, the AO Hessians and the
/// Evaluator results alike, lists the upper triangle row by row: xx, xy, xz,
/// yy, yz, zz.
inline constexpr int ao_hessian_column[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};

} // namespace chem::kumo::detail
