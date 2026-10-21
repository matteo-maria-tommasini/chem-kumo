#include "chem/kumo/local_spatial_persistent_tree.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
using namespace chem::kumo;

static double detailed(const std::vector<GaussianPairEnvelope>& a,
                       const std::vector<GaussianPairEnvelope>& b) {
    if(a.size()==1 && b.size()==1) return local_gaussian_envelope_overlap_bound(a[0],b[0]);
    long double s=0;
    for(const auto& x:a) for(const auto& y:b) s+=local_gaussian_envelope_overlap_bound(x,y);
    if(s>=std::numeric_limits<double>::max()) return INFINITY;
    return std::nextafter(static_cast<double>(s),INFINITY);
}
int main() {
    std::mt19937_64 rng(160316);
    std::uniform_real_distribution<double> coord(-12,12),expo(.2,8),amp(.01,2);
    std::size_t failures=0,better=0,nodes_old=0,nodes_new=0,nodes_budget=0;
    double worst_shortfall=0;
    long double ratios=0;
    double build_seconds=0,old_seconds=0,new_seconds=0,budget_seconds=0;
    for(int trial=0;trial<500;++trial) {
        std::vector<GaussianPairEnvelope> a,b;
        for(int i=0;i<1+trial%13;++i) a.push_back({expo(rng),{coord(rng),coord(rng),coord(rng)},amp(rng)});
        for(int i=0;i<1+(trial*7)%17;++i) b.push_back({expo(rng),{coord(rng),coord(rng),coord(rng)},amp(rng)});
        const auto t0=std::chrono::steady_clock::now();
        SpatialPersistentTree ta(a),tb(b);
        const auto t1=std::chrono::steady_clock::now();
        const auto old=local_spatial_adaptive_bound(a,b,8);
        const auto t2=std::chrono::steady_clock::now();
        const auto now=local_spatial_persistent_bound(ta,tb,8);
        const auto t3=std::chrono::steady_clock::now();
        const double ref=detailed(a,b);
        const double budget=std::isfinite(ref) ? ref*1.1 : INFINITY;
        const auto stopped=local_spatial_persistent_bound(ta,tb,8,budget,true);
        const auto t4=std::chrono::steady_clock::now();
        auto seconds=[](auto x,auto y){return std::chrono::duration<double>(y-x).count();};
        build_seconds+=seconds(t0,t1); old_seconds+=seconds(t1,t2);
        new_seconds+=seconds(t2,t3); budget_seconds+=seconds(t3,t4);
        nodes_old+=old.statistics.node_visits; nodes_new+=now.statistics.node_visits;
        nodes_budget+=stopped.statistics.node_visits;
        if(now.value<old.value) ++better;
        if(now.value<ref || stopped.value<ref || (stopped.within_budget && stopped.value>budget)) {
            ++failures;
            worst_shortfall=std::max(worst_shortfall,(ref-now.value)/std::max(ref,1e-300));
        }
        if(ref>0 && std::isfinite(ref)) ratios+=now.value/ref;
    }
    std::cout<<std::setprecision(12)<<"cases=500 failures="<<failures
        <<" persistent_better="<<better<<" mean_ratio="<<double(ratios/500)
        <<" mean_old_nodes="<<double(nodes_old)/500
        <<" mean_persistent_nodes="<<double(nodes_new)/500
        <<" mean_budget_nodes="<<double(nodes_budget)/500
        <<" build_seconds="<<build_seconds<<" recursive_seconds="<<old_seconds
        <<" persistent_seconds="<<new_seconds<<" budget_seconds="<<budget_seconds
        <<" worst_relative_shortfall="<<worst_shortfall
        <<" formally_certified=0 "<<(failures?"FAIL":"OK")<<'\n';
    return failures?1:0;
}
