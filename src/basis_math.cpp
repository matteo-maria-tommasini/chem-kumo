#include "basis_math.hpp"

#include <cmath>
#include <numbers>

namespace chem::kumo::detail {

double double_factorial_odd(int n) noexcept {
    if (n <= 0) return 1.0;
    double value = 1.0;
    for (int i = n; i >= 1; i -= 2) value *= static_cast<double>(i);
    return value;
}

double cartesian_primitive_normalization(double alpha, const CartesianPower& p) {
    const int l = p.x + p.y + p.z;
    const double radial = std::pow(2.0 * alpha / std::numbers::pi, 0.75);
    const double angular = std::sqrt(std::pow(4.0 * alpha, static_cast<double>(l)) /
        (double_factorial_odd(2 * p.x - 1) *
         double_factorial_odd(2 * p.y - 1) *
         double_factorial_odd(2 * p.z - 1)));
    return radial * angular;
}

} // namespace chem::kumo::detail
