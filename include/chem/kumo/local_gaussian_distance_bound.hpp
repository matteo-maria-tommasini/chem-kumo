#pragma once
#include "chem/kumo/local_gaussian_pair.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace chem::kumo {

// Upper bound on the absolute integral of two primitive Gaussian products.
// All polynomial translations and moment summations use long double to
// reduce roundoff. This is NOT an interval-arithmetic certified bound.
[[nodiscard]] inline double local_gaussian_pair_distance_bound(
    const GaussianPairIntermediate& a, const GaussianPairIntermediate& b) {
    if (!(a.exponent > 0) || !(b.exponent > 0) ||
        !std::isfinite(a.exponent) || !std::isfinite(b.exponent))
        throw std::invalid_argument("invalid pair exponent");
    const long double p = static_cast<long double>(a.exponent) + b.exponent;
    long double r2 = 0;
    std::array<long double, 3> q{};
    for (int d = 0; d < 3; ++d) {
        if (a.degree[d] < 0 || b.degree[d] < 0 ||
            a.degree[d] > local_pair_max_degree ||
            b.degree[d] > local_pair_max_degree)
            throw std::invalid_argument("invalid polynomial degree");
        const long double delta = static_cast<long double>(a.center[d]) - b.center[d];
        r2 += delta * delta;
        q[d] = static_cast<long double>(a.center[d]) +
            (static_cast<long double>(b.exponent) / p) *
            (static_cast<long double>(b.center[d]) - a.center[d]);
    }
    const long double eta = static_cast<long double>(a.exponent) * b.exponent / p;
    long double bound = std::abs(static_cast<long double>(a.prefactor) * b.prefactor)
        * std::exp(-eta * r2);
    const long double pi = std::numbers::pi_v<long double>;
    for (int d = 0; d < 3; ++d) {
        std::array<long double, 2 * local_pair_max_degree + 1> poly{};
        poly[0] = 1;
        int total_degree = 0;
        for (const auto* pair : {&a, &b}) {
            std::array<long double, 2 * local_pair_max_degree + 1> shifted{};
            const long double shift = q[d] - pair->center[d];
            // Horner translation of polynomial in (r - pair.center).
            int degree = 0;
            for (int i = pair->degree[d]; i >= 0; --i) {
                if (i != pair->degree[d]) {
                    for (int j = degree + 1; j >= 0; --j)
                        shifted[j] = (j <= degree ? shift * shifted[j] : 0.L)
                                   + (j > 0 ? shifted[j - 1] : 0.L);
                    ++degree;
                }
                shifted[0] += pair->polynomial[d][i];
            }
            std::array<long double, 2 * local_pair_max_degree + 1> next{};
            for (int i = 0; i <= total_degree; ++i)
                for (int j = 0; j <= pair->degree[d]; ++j)
                    next[i + j] += poly[i] * shifted[j];
            poly = next;
            total_degree += pair->degree[d];
        }
        long double moment = std::sqrt(pi / p);
        long double sum = 0;
        for (int n = 0; n <= total_degree; n += 2) {
            if (n) moment *= (n - 1) / (2 * p);
            sum += std::abs(poly[n]) * moment;
        }
        bound *= sum;
    }
    if (std::isnan(bound)) throw std::runtime_error("NaN Gaussian bound");
    if (bound > std::numeric_limits<double>::max())
        return std::numeric_limits<double>::infinity();
    // Small outward safety margin for the final conversion, not a proof
    // of directed rounding for all intermediate operations.
    const double result = static_cast<double>(bound);
    return std::nextafter(result, std::numeric_limits<double>::infinity());
}

} // namespace chem::kumo
