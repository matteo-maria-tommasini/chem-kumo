#pragma once
#include "chem/kumo/local_spatial_adaptive_bound.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace chem::kumo {

struct SpatialPersistentNode {
    SpatialEnvelopeGroup group{};
    std::size_t begin=0, end=0;
    std::size_t left=0, right=0;
    bool leaf=true;
};

class SpatialPersistentTree {
public:
    explicit SpatialPersistentTree(const std::vector<GaussianPairEnvelope>& leaves)
        : leaves_(leaves), order_(leaves.size()) {
        std::iota(order_.begin(),order_.end(),0);
        if(!leaves.empty()) build(0,leaves.size());
    }
    const SpatialPersistentNode& node(std::size_t i) const { return nodes_.at(i); }
    std::size_t node_count() const { return nodes_.size(); }
    std::size_t leaf_count() const { return leaves_.size(); }
    const GaussianPairEnvelope& leaf(std::size_t pos) const { return leaves_[order_[pos]]; }
private:
    std::vector<GaussianPairEnvelope> leaves_;
    std::vector<std::size_t> order_;
    std::vector<SpatialPersistentNode> nodes_;
    std::size_t build(std::size_t begin,std::size_t end) {
        const std::size_t id=nodes_.size();
        nodes_.emplace_back();
        std::vector<GaussianPairEnvelope> subset;
        subset.reserve(end-begin);
        for(std::size_t k=begin;k<end;++k) subset.push_back(leaves_[order_[k]]);
        nodes_[id].group=make_spatial_envelope_group(subset);
        nodes_[id].begin=begin; nodes_[id].end=end;
        if(end-begin==1) return id;
        std::array<double,3> lo{INFINITY,INFINITY,INFINITY};
        std::array<double,3> hi{-INFINITY,-INFINITY,-INFINITY};
        for(std::size_t k=begin;k<end;++k) for(int d=0;d<3;++d) {
            const auto x=leaves_[order_[k]].center[d];
            lo[d]=std::min(lo[d],x); hi[d]=std::max(hi[d],x);
        }
        int axis=0;
        for(int d=1;d<3;++d) if(hi[d]-lo[d]>hi[axis]-lo[axis]) axis=d;
        std::stable_sort(order_.begin()+begin,order_.begin()+end,
            [&](std::size_t a,std::size_t b){return leaves_[a].center[axis]<leaves_[b].center[axis];});
        const auto mid=begin+(end-begin)/2;
        const auto left=build(begin,mid);
        const auto right=build(mid,end);
        nodes_[id].left=left; nodes_[id].right=right; nodes_[id].leaf=false;
        return id;
    }
};

struct SpatialPersistentQuery {
    double value=0;
    SpatialAdaptiveStatistics statistics{};
    bool within_budget=false;
    bool formally_certified=false;
};

namespace spatial_persistent_detail {
inline double evaluate(const SpatialPersistentTree& a,std::size_t ia,
    const SpatialPersistentTree& b,std::size_t ib,std::size_t depth,
    std::size_t max_depth,double budget,bool stop_on_budget,
    SpatialAdaptiveStatistics& stats) {
    ++stats.node_visits;
    const auto& na=a.node(ia);
    const auto& nb=b.node(ib);
    const double parent=local_spatial_aggregate_bound_guarded(na.group,nb.group).value;
    if(stop_on_budget && parent<=budget) return parent;
    if(depth>=max_depth || (na.leaf && nb.leaf)) {
        if(na.leaf && nb.leaf) ++stats.leaf_pairs;
        return parent;
    }
    const bool split_a=!na.leaf && (nb.leaf || na.group.radius>nb.group.radius ||
        (na.group.radius==nb.group.radius && na.end-na.begin>=nb.end-nb.begin));
    ++stats.splits;
    const double first=split_a ? evaluate(a,na.left,b,ib,depth+1,max_depth,budget,stop_on_budget,stats)
                                 : evaluate(a,ia,b,nb.left,depth+1,max_depth,budget,stop_on_budget,stats);
    // Once the first child alone exceeds the budget, no need to evaluate the
    // second child for a negative decision. Return the parent (a valid bound).
    if(stop_on_budget && first>budget) return parent;
    const double second=split_a ? evaluate(a,na.right,b,ib,depth+1,max_depth,
                                            std::max(0.0,budget-first),stop_on_budget,stats)
                                  : evaluate(a,ia,b,nb.right,depth+1,max_depth,
                                             std::max(0.0,budget-first),stop_on_budget,stats);
    return std::min(parent,spatial_adaptive_detail::add_upper(first,second));
}
} // namespace spatial_persistent_detail

[[nodiscard]] inline SpatialPersistentQuery local_spatial_persistent_bound(
    const SpatialPersistentTree& a,const SpatialPersistentTree& b,
    std::size_t max_depth=8,double budget=0,bool stop_on_budget=false) {
    if(!(budget>=0) || std::isnan(budget)) throw std::invalid_argument("invalid spatial budget");
    SpatialPersistentQuery out;
    if(a.leaf_count()==0 || b.leaf_count()==0) {out.within_budget=true;return out;}
    out.value=spatial_persistent_detail::evaluate(a,0,b,0,0,max_depth,budget,stop_on_budget,out.statistics);
    out.within_budget=out.value<=budget;
    return out;
}
} // namespace chem::kumo
