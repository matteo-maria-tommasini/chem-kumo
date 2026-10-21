#pragma once
#include "chem/kumo/local_spatial_aggregate_bound.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace chem::kumo {

struct SpatialAdaptiveStatistics {
    std::size_t node_visits=0;
    std::size_t splits=0;
    std::size_t leaf_pairs=0;
};

struct SpatialAdaptiveResult {
    double value=0;
    SpatialAdaptiveStatistics statistics{};
    bool formally_certified=false;
};

namespace spatial_adaptive_detail {
inline double add_upper(double x,double y) {
    if(std::isinf(x)||std::isinf(y)) return std::numeric_limits<double>::infinity();
    return std::nextafter(x+y,std::numeric_limits<double>::infinity());
}
inline std::pair<std::vector<GaussianPairEnvelope>,std::vector<GaussianPairEnvelope>> split(
    const std::vector<GaussianPairEnvelope>& in) {
    std::array<double,3> lo{INFINITY,INFINITY,INFINITY};
    std::array<double,3> hi{-INFINITY,-INFINITY,-INFINITY};
    for(const auto& e:in) for(int d=0;d<3;++d) {
        lo[d]=std::min(lo[d],e.center[d]); hi[d]=std::max(hi[d],e.center[d]);
    }
    int axis=0;
    for(int d=1;d<3;++d) if(hi[d]-lo[d]>hi[axis]-lo[axis]) axis=d;
    auto sorted=in;
    std::stable_sort(sorted.begin(),sorted.end(),[axis](const auto& a,const auto& b){return a.center[axis]<b.center[axis];});
    const auto mid=sorted.size()/2;
    return {{sorted.begin(),sorted.begin()+mid},{sorted.begin()+mid,sorted.end()}};
}
inline double evaluate(const std::vector<GaussianPairEnvelope>& a,
    const std::vector<GaussianPairEnvelope>& b,std::size_t depth,
    std::size_t max_depth,SpatialAdaptiveStatistics& stats) {
    ++stats.node_visits;
    const auto ga=make_spatial_envelope_group(a);
    const auto gb=make_spatial_envelope_group(b);
    const double parent=local_spatial_aggregate_bound_guarded(ga,gb).value;
    if(a.empty()||b.empty()) return 0;
    if(depth>=max_depth || (a.size()==1 && b.size()==1)) {
        if(a.size()==1 && b.size()==1) ++stats.leaf_pairs;
        return parent;
    }
    // Split the group with greater spatial extent; if tied, split the larger group.
    const bool split_a=(a.size()>1) && (b.size()==1 || ga.radius>gb.radius ||
        (ga.radius==gb.radius && a.size()>=b.size()));
    if(split_a) {
        auto [left,right]=split(a); ++stats.splits;
        return std::min(parent,add_upper(evaluate(left,b,depth+1,max_depth,stats),
                                          evaluate(right,b,depth+1,max_depth,stats)));
    }
    auto [left,right]=split(b); ++stats.splits;
    return std::min(parent,add_upper(evaluate(a,left,depth+1,max_depth,stats),
                                      evaluate(a,right,depth+1,max_depth,stats)));
}
} // namespace spatial_adaptive_detail

[[nodiscard]] inline SpatialAdaptiveResult local_spatial_adaptive_bound(
    const std::vector<GaussianPairEnvelope>& a,
    const std::vector<GaussianPairEnvelope>& b,std::size_t max_depth=8) {
    SpatialAdaptiveResult result;
    result.value=spatial_adaptive_detail::evaluate(a,b,0,max_depth,result.statistics);
    return result;
}
} // namespace chem::kumo
