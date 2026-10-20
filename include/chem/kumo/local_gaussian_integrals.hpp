#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>
#include <utility>
#include <cstddef>

namespace chem::kumo {

/// Unnormalized Cartesian primitive; scale includes any primitive normalization
/// and contraction coefficient. Coordinates are in bohr; exponents in bohr^-2.
struct LocalGaussianPrimitive {
    std::array<double, 3> center{};
    std::array<int, 3> powers{};
    double exponent = 1.0;
    double scale = 1.0;
};

namespace detail {
inline std::vector<double> multiply_shifted_power(
    const std::vector<double>& polynomial, double shift, int power) {
    std::vector<double> result = polynomial;
    for (int j = 0; j < power; ++j) {
        std::vector<double> next(result.size() + 1, 0.0);
        for (std::size_t k = 0; k < result.size(); ++k) {
            next[k] += shift * result[k];
            next[k + 1] += result[k];
        }
        result = std::move(next);
    }
    return result;
}

inline double gaussian_polynomial_moment(
    const std::array<LocalGaussianPrimitive, 4>& g,
    int axis, double center, double exponent_sum) {
    std::vector<double> polynomial{1.0};
    for (const auto& primitive : g)
        polynomial = multiply_shifted_power(polynomial,
            center - primitive.center[static_cast<std::size_t>(axis)],
            primitive.powers[static_cast<std::size_t>(axis)]);
    const double base = std::sqrt(std::numbers::pi / exponent_sum);
    double moment = base;
    double sum = 0.0;
    for (std::size_t n = 0; n < polynomial.size(); n += 2) {
        if (n != 0) moment *= static_cast<double>(n - 1) / (2.0 * exponent_sum);
        sum += polynomial[n] * moment;
    }
    return sum;
}
} // namespace detail

/// Exact one-coordinate-space integral of four Cartesian Gaussian primitives.
/// No Coulomb operator, no four-index tensor, no numerical quadrature.
[[nodiscard]] inline double local_four_gaussian_integral(
    const std::array<LocalGaussianPrimitive, 4>& g) {
    double p = 0.0;
    std::array<double, 3> q{};
    double scale = 1.0;
    for (const auto& primitive : g) {
        if (!(primitive.exponent > 0.0) || !std::isfinite(primitive.exponent)
            || !std::isfinite(primitive.scale))
            throw std::invalid_argument("invalid Gaussian exponent or scale");
        p += primitive.exponent;
        scale *= primitive.scale;
        for (int d = 0; d < 3; ++d) {
            if (primitive.powers[static_cast<std::size_t>(d)] < 0 ||
                !std::isfinite(primitive.center[static_cast<std::size_t>(d)]))
                throw std::invalid_argument("invalid Gaussian center or power");
            q[static_cast<std::size_t>(d)] +=
                primitive.exponent * primitive.center[static_cast<std::size_t>(d)];
        }
    }
    if (!std::isfinite(p)) throw std::invalid_argument("invalid exponent sum");
    for (auto& x : q) x /= p;
    // Pairwise form avoids cancellation between large absolute center coordinates.
    double separation = 0.0;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = i + 1; j < 4; ++j) {
            double r2 = 0.0;
            for (int d = 0; d < 3; ++d) {
                const double delta = g[i].center[static_cast<std::size_t>(d)] -
                                     g[j].center[static_cast<std::size_t>(d)];
                r2 += delta * delta;
            }
            separation += (g[i].exponent / p) * g[j].exponent * r2;
        }
    double value = scale * std::exp(-separation);
    for (int d = 0; d < 3; ++d)
        value *= detail::gaussian_polynomial_moment(g, d,
            q[static_cast<std::size_t>(d)], p);
    return value;
}
} // namespace chem::kumo
