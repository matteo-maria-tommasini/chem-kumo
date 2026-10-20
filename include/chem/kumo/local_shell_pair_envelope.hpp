#pragma once

#include "chem/kumo/local_gaussian_pair_envelope.hpp"
#include "chem/kumo/local_shell_integrals.hpp"
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

// An envelope for every AO-component product in a shell pair. All primitive
// coefficients and Cartesian-to-pure transformation factors are already
// included in LocalGaussianPrimitive::scale.
struct LocalShellPairEnvelope {
    std::size_t na = 0, nb = 0;
    std::vector<std::vector<GaussianPairEnvelope>> components;
    [[nodiscard]] const std::vector<GaussianPairEnvelope>& at(
        std::size_t i, std::size_t j) const {
        return components.at(i * nb + j);
    }
    [[nodiscard]] std::size_t envelope_count() const {
        std::size_t n = 0;
        for (const auto& v : components) n += v.size();
        return n;
    }
};

[[nodiscard]] inline LocalShellPairEnvelope prepare_local_shell_pair_envelope(
    const LocalShellComponents& a, const LocalShellComponents& b,
    double gamma = 0.7) {
    LocalShellPairEnvelope result;
    result.na = a.size();
    result.nb = b.size();
    result.components.reserve(result.na * result.nb);
    for (const auto& ca : a)
        for (const auto& cb : b) {
            std::vector<GaussianPairEnvelope> env;
            env.reserve(ca.size() * cb.size());
            for (const auto& pa : ca)
                for (const auto& pb : cb)
                    env.push_back(make_gaussian_pair_envelope(
                        make_gaussian_pair(pa, pb), gamma));
            result.components.push_back(std::move(env));
        }
    return result;
}

// Componentwise bound for the integral of two contracted AO products.
// Triangle inequality applies across both primitive expansions.
[[nodiscard]] inline double local_shell_pair_envelope_bound(
    const LocalShellPairEnvelope& ab, std::size_t i, std::size_t j,
    const LocalShellPairEnvelope& cd, std::size_t k, std::size_t l) {
    long double bound = 0.0L;
    for (const auto& a : ab.at(i, j))
        for (const auto& b : cd.at(k, l))
            bound += local_gaussian_envelope_overlap_bound(a, b);
    if (bound >= std::numeric_limits<double>::max())
        return std::numeric_limits<double>::infinity();
    return std::nextafter(static_cast<double>(bound),
                          std::numeric_limits<double>::infinity());
}

// Maximum elementwise bound over a complete shell quartet.
[[nodiscard]] inline double local_shell_quartet_envelope_bound(
    const LocalShellPairEnvelope& ab, const LocalShellPairEnvelope& cd) {
    double bound = 0.0;
    for (std::size_t i = 0; i < ab.na; ++i)
        for (std::size_t j = 0; j < ab.nb; ++j)
            for (std::size_t k = 0; k < cd.na; ++k)
                for (std::size_t l = 0; l < cd.nb; ++l)
                    bound = std::max(bound,
                        local_shell_pair_envelope_bound(ab, i, j, cd, k, l));
    return bound;
}

} // namespace chem::kumo
