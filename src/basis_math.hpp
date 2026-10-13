#pragma once

#include "chem/kumo/basis_evaluator.hpp"

namespace chem::kumo::detail {

double double_factorial_odd(int n) noexcept;
double cartesian_primitive_normalization(double alpha, const CartesianPower& power);

} // namespace chem::kumo::detail
