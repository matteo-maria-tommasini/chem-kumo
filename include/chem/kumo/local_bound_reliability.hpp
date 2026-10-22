#pragma once
#include <cmath>
#include <limits>

namespace chem::kumo {

enum class BoundReliability { Certified, Engineering, Unusable };
enum class BoundPolicy { Rigorous, Experimental };

struct ReliableBound {
    double value = std::numeric_limits<double>::infinity();
    BoundReliability reliability = BoundReliability::Unusable;
};

// This gate does not certify any bound. Certification must originate
// from a separately validated, outward-rounded bound implementation.
[[nodiscard]] inline bool local_bound_allows_skip(
    ReliableBound bound, double budget,
    BoundPolicy policy = BoundPolicy::Rigorous) noexcept {
    if (!std::isfinite(budget) || budget <= 0.0 ||
        !std::isfinite(bound.value) || bound.value < 0.0 ||
        bound.reliability == BoundReliability::Unusable)
        return false;
    if (policy == BoundPolicy::Rigorous &&
        bound.reliability != BoundReliability::Certified)
        return false;
    // A computed zero is not evidence of exact zero; reject it even
    // when the caller labels the bound certified.
    if (bound.value == 0.0) return false;
    return bound.value <= budget;
}

} // namespace chem::kumo
