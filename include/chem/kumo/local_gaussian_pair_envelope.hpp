#pragma once
#include "chem/kumo/local_gaussian_pair.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace chem::kumo {

struct GaussianPairEnvelope {
    double exponent = 0.0;
    std::array<double, 3> center{};
    double amplitude = 0.0;
};

// Pointwise envelope: |pair(r)| <= amplitude * exp(-exponent*|r-center|^2).
// The polynomial is bounded after reserving (1-gamma) of the exponent.
// Floating-point evaluation is not a formally certified outward bound.
[[nodiscard]] inline GaussianPairEnvelope make_gaussian_pair_envelope(
    const GaussianPairIntermediate& pair, double gamma = 0.5) {
    if (!(gamma > 0.0 && gamma < 1.0) || !std::isfinite(gamma) ||
        !(pair.exponent > 0.0) || !std::isfinite(pair.exponent))
        throw std::invalid_argument("invalid Gaussian envelope parameters");
    const long double beta = (1.0L - gamma) * pair.exponent;
    long double amplitude = std::abs(static_cast<long double>(pair.prefactor));
    for (int d = 0; d < 3; ++d) {
        if (pair.degree[d] < 0 || pair.degree[d] > local_pair_max_degree)
            throw std::invalid_argument("invalid Gaussian pair polynomial degree");
        long double axis = 0.0L;
        for (int n = 0; n <= pair.degree[d]; ++n) {
            const long double c = std::abs(static_cast<long double>(pair.polynomial[d][n]));
            const long double supremum = n == 0 ? 1.0L :
                std::pow(static_cast<long double>(n) / (2.0L * std::numbers::e_v<long double> * beta),
                         static_cast<long double>(n) / 2.0L);
            axis += c * supremum;
        }
        amplitude *= axis;
    }
    if (std::isnan(amplitude)) throw std::runtime_error("NaN envelope amplitude");
    const double a = amplitude >= std::numeric_limits<double>::max() ?
        std::numeric_limits<double>::infinity() :
        std::nextafter(static_cast<double>(amplitude), std::numeric_limits<double>::infinity());
    return {gamma * pair.exponent, pair.center, a};
}

// Analytic overlap of the two positive Gaussian envelopes.
[[nodiscard]] inline double local_gaussian_envelope_overlap_bound(
    const GaussianPairEnvelope& a, const GaussianPairEnvelope& b) {
    if (!(a.exponent > 0.0 && b.exponent > 0.0) ||
        !(a.amplitude >= 0.0 && b.amplitude >= 0.0))
        throw std::invalid_argument("invalid Gaussian envelope");
    const long double p = static_cast<long double>(a.exponent) + b.exponent;
    const long double eta = static_cast<long double>(a.exponent) * b.exponent / p;
    long double distance2 = 0;
    for (int d = 0; d < 3; ++d) {
        const long double delta = static_cast<long double>(a.center[d]) - b.center[d];
        distance2 += delta * delta;
    }
    const long double result = static_cast<long double>(a.amplitude) * b.amplitude *
        std::pow(std::numbers::pi_v<long double> / p, 1.5L) * std::exp(-eta * distance2);
    if (std::isnan(result)) throw std::runtime_error("NaN envelope overlap");
    if (result >= std::numeric_limits<double>::max())
        return std::numeric_limits<double>::infinity();
    return std::nextafter(static_cast<double>(result), std::numeric_limits<double>::infinity());
}

} // namespace chem::kumo
