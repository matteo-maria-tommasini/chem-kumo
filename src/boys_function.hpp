#pragma once
#include <span>

namespace chem::kumo::detail {
/// Largest order m accepted by boys_function(): enough for Coulomb integrals
/// of shell pairs with total angular momentum up to 32.
inline constexpr int boys_max_order = 32;

/// Boys function F_m(T) = integral_0^1 t^(2m) exp(-T t^2) dt for
/// m = 0, ..., values.size() - 1, written to `values`.
///
/// For T < 40, F_M (M the highest order requested) comes from a seven-term
/// (sixth-degree) Taylor expansion about the nearest point of a table with spacing 0.05, and
/// the lower orders from the downward recursion
///   F_(m-1) = (2 T F_m + exp(-T)) / (2m - 1),
/// which is stable. For T >= 40, F_0 = sqrt(pi / T) / 2 (erf(sqrt(T)) differs
/// from 1 by less than 1e-18) and the upward recursion
///   F_(m+1) = ((2m + 1) F_m - exp(-T)) / (2 T)
/// is stable because 2m + 1 < 2T for every supported order. The relative
/// error is about 1e-15 (see tests/boys_function_test.cpp).
///
/// Throws std::invalid_argument for an empty span, more than
/// boys_max_order + 1 values, or T negative or not finite.
void boys_function(double T, std::span<double> values);
} // namespace chem::kumo::detail
