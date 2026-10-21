#pragma once
#include "chem/kumo/local_gaussian_pair_envelope.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

// Engineering-level guard, NOT a formally outward-rounded interval bound.
// In particular the input envelopes and libm transcendental functions are not certified.
struct SpatialEnvelopeGroup {
    std::array<long double,3> center{};
    long double radius=0;
    long double minimum_exponent=std::numeric_limits<long double>::infinity();
    long double amplitude_sum=0;
    std::size_t size=0;
    GaussianPairEnvelope singleton{}; // Valid only when size == 1.
};

struct SpatialAggregateBoundResult {
    double value=0;
    long double unguarded=0;
    long double guard=0;
    bool formally_certified=false;
};

[[nodiscard]] inline SpatialEnvelopeGroup make_spatial_envelope_group(
    const std::vector<GaussianPairEnvelope>& leaves) {
    SpatialEnvelopeGroup g;
    g.size=leaves.size();
    if(leaves.empty()) return g;
    if(leaves.size()==1) g.singleton=leaves.front();
    for(const auto& e:leaves) {
        if(!(e.exponent>0) || !std::isfinite(e.exponent) ||
           !(e.amplitude>=0) || !std::isfinite(e.amplitude))
            throw std::invalid_argument("invalid spatial envelope");
        for(int d=0;d<3;++d) {
            if(!std::isfinite(e.center[d])) throw std::invalid_argument("invalid envelope center");
            g.center[d]+=static_cast<long double>(e.center[d])/leaves.size();
        }
        g.minimum_exponent=std::min(g.minimum_exponent,static_cast<long double>(e.exponent));
        g.amplitude_sum+=static_cast<long double>(e.amplitude);
    }
    for(const auto& e:leaves) {
        long double r2=0;
        for(int d=0;d<3;++d) {
            const auto dx=static_cast<long double>(e.center[d])-g.center[d];
            r2+=dx*dx;
        }
        g.radius=std::max(g.radius,std::sqrt(r2));
    }
    return g;
}

[[nodiscard]] inline SpatialAggregateBoundResult local_spatial_aggregate_bound_guarded(
    const SpatialEnvelopeGroup& a,const SpatialEnvelopeGroup& b) {
    if(a.size==0 || b.size==0) return {};
    if(!(a.minimum_exponent>0 && b.minimum_exponent>0) ||
       !(a.amplitude_sum>=0 && b.amplitude_sum>=0) ||
       !(a.radius>=0 && b.radius>=0))
        throw std::invalid_argument("invalid spatial envelope group");
    // Exact API-consistent singleton fast path; no duplicate transcendental evaluation.
    if(a.size==1 && b.size==1) {
        const double v=local_gaussian_envelope_overlap_bound(a.singleton,b.singleton);
        return {v,static_cast<long double>(v),0.0L,false};
    }
    long double d2=0;
    for(int k=0;k<3;++k) {
        const auto delta=a.center[k]-b.center[k];
        d2+=delta*delta;
    }
    const auto separation=std::max(0.0L,std::sqrt(d2)-a.radius-b.radius);
    const auto p=a.minimum_exponent+b.minimum_exponent;
    const auto eta=a.minimum_exponent*b.minimum_exponent/p;
    const auto pi=std::acos(-1.0L);
    const auto raw=a.amplitude_sum*b.amplitude_sum*std::pow(pi/p,1.5L)*
        std::exp(-eta*separation*separation);
    if(std::isnan(raw) || raw<0) throw std::runtime_error("invalid aggregate bound");
    // Deliberately conservative engineering allowance for mixed double/long-double
    // computations and libm rounding. This is NOT a proof of a floating-point bound.
    constexpr long double relative_guard=128.0L*std::numeric_limits<double>::epsilon();
    const auto guard=raw*relative_guard;
    const auto padded=raw+guard;
    const double value=padded>=std::numeric_limits<double>::max() ?
        std::numeric_limits<double>::infinity() :
        std::nextafter(static_cast<double>(padded),std::numeric_limits<double>::infinity());
    return {value,raw,guard,false};
}
} // namespace chem::kumo
