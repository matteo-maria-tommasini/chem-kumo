#pragma once
// Fail-closed adapter for the optional experimental MPFR s-quartet kernel.
// Does not confer mathematical certification on the underlying bound.
#include "chem/kumo/local_bound_reliability.hpp"
#include "chem/kumo/local_mpfr_s_quartet_bound.hpp"
#include <cmath>
#include <limits>

namespace chem::kumo::certified {
enum class SQuartetBoundStatus {
    FiniteEngineering,
    ExactZeroCoefficient,
    InvalidInput,
    ConservativeFallback
};
struct SQuartetBoundResult {
    chem::kumo::ReliableBound bound{};
    SQuartetBoundStatus status = SQuartetBoundStatus::ConservativeFallback;
};

[[nodiscard]] inline SQuartetBoundResult safe_s_primitive_quartet_bound(
    const SPrimitive& a, const SPrimitive& b,
    const SPrimitive& c, const SPrimitive& d) noexcept {
    const SPrimitive* p[] = {&a,&b,&c,&d};
    for (const auto* v : p) {
        if (!(v->exponent > 0.0) || !std::isfinite(v->exponent) ||
            !std::isfinite(v->coefficient))
            return {{}, SQuartetBoundStatus::InvalidInput};
        for (double x : v->center)
            if (!std::isfinite(x))
                return {{}, SQuartetBoundStatus::InvalidInput};
    }
    for (const auto* v : p)
        if (v->coefficient == 0.0)
            return {{0.0, chem::kumo::BoundReliability::Engineering},
                    SQuartetBoundStatus::ExactZeroCoefficient};
    try {
        const double upper = s_primitive_quartet_overlap_upper(a,b,c,d);
        if (!std::isfinite(upper) || !(upper > 0.0))
            return {{}, SQuartetBoundStatus::ConservativeFallback};
        return {{upper, chem::kumo::BoundReliability::Engineering},
                SQuartetBoundStatus::FiniteEngineering};
    } catch (...) {
        return {{}, SQuartetBoundStatus::ConservativeFallback};
    }
}
} // namespace chem::kumo::certified
