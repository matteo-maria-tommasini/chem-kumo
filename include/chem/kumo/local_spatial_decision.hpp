#pragma once
#include "chem/kumo/local_spatial_persistent_tree.hpp"
#include <cmath>
#include <cstddef>
#include <limits>

namespace chem::kumo {

// CertifiedSkip means accepted under the *engineering* floating-point bound.
// It does not imply formally verified outward rounding.
enum class SpatialDecision { CertifiedSkip, KeepOrRefine };

struct SpatialDecisionResult {
    SpatialDecision decision=SpatialDecision::KeepOrRefine;
    double upper_bound=std::numeric_limits<double>::infinity();
    std::size_t node_visits=0;
    bool formally_certified=false;
};

[[nodiscard]] inline SpatialDecisionResult local_spatial_persistent_decision(
    const SpatialPersistentTree& a, const SpatialPersistentTree& b,
    std::size_t max_depth, double budget) {
    SpatialDecisionResult out;
    // Reject invalid budgets conservatively; do not call the underlying query.
    if (!std::isfinite(budget) || budget < 0.0) return out;
    // An empty set has identically zero overlap.
    if (a.leaf_count()==0 || b.leaf_count()==0) {
        out.upper_bound=0.0;
        out.decision=SpatialDecision::CertifiedSkip;
        return out;
    }
    // A zero budget must not permit skipping a positive (or underflowed) bound.
    // In particular, do not trust a floating-point zero from exp underflow.
    if (budget==0.0) return out;
    const auto query=local_spatial_persistent_bound(a,b,max_depth,budget,true);
    out.node_visits=query.statistics.node_visits;
    if (!std::isfinite(query.value) || query.value < 0.0) return out;
    out.upper_bound=query.value;
    if (query.within_budget && query.value<=budget)
        out.decision=SpatialDecision::CertifiedSkip;
    return out;
}

} // namespace chem::kumo
