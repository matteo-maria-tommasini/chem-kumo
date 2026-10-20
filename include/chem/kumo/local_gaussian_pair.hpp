#pragma once

#include "chem/kumo/local_gaussian_integrals.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace chem::kumo {

// Pair polynomial is expressed about the Gaussian product center.
// Fixed capacity supports total Cartesian degree <= 12 per axis.
inline constexpr int local_pair_max_degree = 12;

struct GaussianPairIntermediate {
    double exponent = 0.0;
    std::array<double, 3> center{};
    double prefactor = 0.0;
    std::array<int, 3> degree{};
    std::array<std::array<double, local_pair_max_degree + 1>, 3> polynomial{};
};

[[nodiscard]] inline GaussianPairIntermediate make_gaussian_pair(
    const LocalGaussianPrimitive& a, const LocalGaussianPrimitive& b) {
    auto valid = [](const LocalGaussianPrimitive& g) {
        if (!(g.exponent > 0.0) || !std::isfinite(g.exponent) ||
            !std::isfinite(g.scale)) return false;
        for (int d = 0; d < 3; ++d)
            if (g.powers[d] < 0 || !std::isfinite(g.center[d])) return false;
        return true;
    };
    if (!valid(a) || !valid(b))
        throw std::invalid_argument("invalid Gaussian primitive");
    GaussianPairIntermediate pair;
    pair.exponent = a.exponent + b.exponent;
    if (!std::isfinite(pair.exponent))
        throw std::invalid_argument("invalid pair exponent");
    double distance2 = 0.0;
    for (int d = 0; d < 3; ++d) {
        const int degree = a.powers[d] + b.powers[d];
        if (degree > local_pair_max_degree || degree < 0)
            throw std::invalid_argument("Gaussian pair angular degree exceeds capacity");
        pair.degree[d] = degree;
        pair.center[d] = a.center[d] +
            (b.exponent / pair.exponent) * (b.center[d] - a.center[d]);
        const double delta = a.center[d] - b.center[d];
        distance2 += delta * delta;
        auto& coefficients = pair.polynomial[d];
        coefficients[0] = 1.0;
        int current_degree = 0;
        for (int primitive = 0; primitive < 2; ++primitive) {
            const int power = primitive == 0 ? a.powers[d] : b.powers[d];
            const double shift = pair.center[d] -
                (primitive == 0 ? a.center[d] : b.center[d]);
            for (int j = 0; j < power; ++j) {
                for (int k = current_degree + 1; k >= 0; --k)
                    coefficients[k] = (k <= current_degree ? shift * coefficients[k] : 0.0)
                        + (k > 0 ? coefficients[k - 1] : 0.0);
                ++current_degree;
            }
        }
    }
    pair.prefactor = a.scale * b.scale *
        std::exp(-(a.exponent / pair.exponent) * b.exponent * distance2);
    return pair;
}

[[nodiscard]] inline double local_gaussian_pair_integral(
    const GaussianPairIntermediate& a, const GaussianPairIntermediate& b) {
    const double exponent = a.exponent + b.exponent;
    double separation = 0.0;
    double value = a.prefactor * b.prefactor;
    for (int d = 0; d < 3; ++d) {
        const double delta = a.center[d] - b.center[d];
        separation += (a.exponent / exponent) * b.exponent * delta * delta;
    }
    value *= std::exp(-separation);
    const double base = std::sqrt(std::numbers::pi / exponent);
    for (int d = 0; d < 3; ++d) {
        const double product_center = a.center[d] +
            (b.exponent / exponent) * (b.center[d] - a.center[d]);
        std::array<double, 2 * local_pair_max_degree + 1> coefficients{};
        coefficients[0] = 1.0;
        int total_degree = 0;
        for (const auto* pair : {&a, &b}) {
            // Horner translation: P(t + shift), not P(t)*(t+shift)^degree.
            std::array<double, 2 * local_pair_max_degree + 1> next{};
            int shifted_degree = 0;
            const double shift = product_center - pair->center[d];
            for (int i = pair->degree[d]; i >= 0; --i) {
                if (i != pair->degree[d]) {
                    for (int k = shifted_degree + 1; k >= 0; --k)
                        next[k] = (k <= shifted_degree ? shift * next[k] : 0.0)
                            + (k > 0 ? next[k - 1] : 0.0);
                    ++shifted_degree;
                }
                next[0] += pair->polynomial[d][i];
            }
            // Shifted polynomial for this pair is now 'next'. Multiply it into
            // the accumulated polynomial (which already includes earlier pairs).
            std::array<double, 2 * local_pair_max_degree + 1> combined{};
            for (int i = 0; i <= total_degree; ++i)
                for (int j = 0; j <= pair->degree[d]; ++j)
                    combined[i + j] += coefficients[i] * next[j];
            coefficients = combined;
            total_degree += pair->degree[d];
        }
        double moment = base;
        double sum = 0.0;
        for (int n = 0; n <= total_degree; n += 2) {
            if (n != 0) moment *= (n - 1) / (2.0 * exponent);
            sum += coefficients[n] * moment;
        }
        value *= sum;
    }
    return value;
}

} // namespace chem::kumo
