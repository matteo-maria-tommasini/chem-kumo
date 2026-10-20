#pragma once

#include "chem/kumo/local_gaussian_pair.hpp"
#include "chem/kumo/local_shell_integrals.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

// One prepared shell pair. Entries are ordered as AO-component i, j.
// Only primitive-pair intermediates are cached, never four-index integrals.
struct LocalPreparedShellPair {
    std::size_t na = 0;
    std::size_t nb = 0;
    std::vector<std::vector<GaussianPairIntermediate>> components;

    [[nodiscard]] const std::vector<GaussianPairIntermediate>& at(
        std::size_t i, std::size_t j) const {
        return components.at(i * nb + j);
    }
    [[nodiscard]] std::size_t pair_count() const {
        std::size_t count = 0;
        for (const auto& entry : components) count += entry.size();
        return count;
    }
};

[[nodiscard]] inline LocalPreparedShellPair prepare_local_shell_pair(
    const LocalShellComponents& a, const LocalShellComponents& b) {
    LocalPreparedShellPair prepared;
    prepared.na = a.size();
    prepared.nb = b.size();
    prepared.components.reserve(prepared.na * prepared.nb);
    for (const auto& ca : a)
        for (const auto& cb : b) {
            std::vector<GaussianPairIntermediate> pairs;
            pairs.reserve(ca.size() * cb.size());
            for (const auto& pa : ca)
                for (const auto& pb : cb)
                    pairs.push_back(make_gaussian_pair(pa, pb));
            prepared.components.push_back(std::move(pairs));
        }
    return prepared;
}

// A quartet of two precomputed shell pairs, with the same component order
// as local_shell_quartet({A,B,C,D}).
[[nodiscard]] inline LocalShellQuartet local_shell_quartet_paired(
    const LocalPreparedShellPair& ab, const LocalPreparedShellPair& cd) {
    LocalShellQuartet result;
    result.dimensions = {ab.na, ab.nb, cd.na, cd.nb};
    result.values.reserve(ab.na * ab.nb * cd.na * cd.nb);
    for (std::size_t i = 0; i < ab.na; ++i)
        for (std::size_t j = 0; j < ab.nb; ++j)
            for (std::size_t k = 0; k < cd.na; ++k)
                for (std::size_t l = 0; l < cd.nb; ++l) {
                    double value = 0.0;
                    for (const auto& p : ab.at(i, j))
                        for (const auto& q : cd.at(k, l))
                            value += local_gaussian_pair_integral(p, q);
                    result.values.push_back(value);
                }
    return result;
}

} // namespace chem::kumo
