#pragma once
#include "chem/kumo/local_spatial_persistent_tree.hpp"
#include "chem/kumo/local_bound_reliability.hpp"
#include <cmath>
#include <cstddef>
#include <limits>

namespace chem::kumo {

// CertifiedSkip means accepted under the *engineering* floating-point bound.
// It does not imply formally verified outward rounding.
enum class SpatialDecision { Skip, KeepOrRefine, CertifiedSkip = Skip };

struct SpatialDecisionResult {
    SpatialDecision decision=SpatialDecision::KeepOrRefine;
    double upper_bound=std::numeric_limits<double>::infinity();
    std::size_t node_visits=0;
    bool formally_certified=false;
};

[[nodiscard]] inline SpatialDecisionResult local_spatial_persistent_decision(
    const SpatialPersistentTree& a, const SpatialPersistentTree& b,
    std::size_t max_depth, double budget, BoundPolicy policy) {
    SpatialDecisionResult out;
    // Reject invalid budgets conservatively; do not call the underlying query.
    if (!std::isfinite(budget) || budget < 0.0) return out;
    // An empty set has identically zero overlap.
    if (a.leaf_count()==0 || b.leaf_count()==0) {
        out.upper_bound=0.0;
        // An empty product is mathematically zero, independently of
        // floating-point envelope evaluation.
        out.formally_certified=true;
        out.decision=SpatialDecision::Skip;
        return out;
    }
    // A zero budget must not permit skipping a positive (or underflowed) bound.
    // In particular, do not trust a floating-point zero from exp underflow.
    if (budget==0.0) return out;
    const auto query=local_spatial_persistent_bound(a,b,max_depth,budget,true);
    out.node_visits=query.statistics.node_visits;
    if (!std::isfinite(query.value) || query.value < 0.0) return out;
    out.upper_bound=query.value;
    // The persistent spatial tree currently supplies engineering bounds only.
    // Never promote them to Certified based on numeric agreement alone.
    const ReliableBound candidate{query.value,BoundReliability::Engineering};
    if (query.within_budget &&
        local_bound_allows_skip(candidate,budget,policy))
        out.decision=SpatialDecision::Skip;
    return out;
}

// Compatibility entry point: retains the historical experimental semantics.
// New safety-sensitive callers must pass BoundPolicy::Rigorous explicitly.
[[deprecated("Pass BoundPolicy explicitly; legacy uses Experimental")]]
[[nodiscard]] inline SpatialDecisionResult local_spatial_persistent_decision(
    const SpatialPersistentTree& a, const SpatialPersistentTree& b,
    std::size_t max_depth, double budget) {
    return local_spatial_persistent_decision(
        a,b,max_depth,budget,BoundPolicy::Experimental);
}

} // namespace chem::kumo
