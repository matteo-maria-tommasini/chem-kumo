#include "chem/kumo/local_spatial_persistent_tree.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
using Clock=std::chrono::steady_clock;

static std::vector<GaussianPairEnvelope> sample(std::mt19937_64& rng,std::size_t n) {
    std::uniform_real_distribution<double> xyz(-12.,12.), alpha(.2,8.), amp(.01,2.);
    std::vector<GaussianPairEnvelope> v;
    v.reserve(n);
    for(std::size_t i=0;i<n;++i) v.push_back({alpha(rng),{xyz(rng),xyz(rng),xyz(rng)},amp(rng)});
    return v;
}
static double seconds(Clock::time_point a,Clock::time_point b) {
    return std::chrono::duration<double>(b-a).count();
}
static bool close(double a,double b) {
    if(a==b) return true;
    return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=1e-12*std::max({1.,std::abs(a),std::abs(b)});
}
int main() {
    try {
        constexpr std::array<std::size_t,6> depths{2,4,6,8,10,12};
        constexpr std::array<double,8> fractions{0.,.01,.1,.5,.9,.99,1.,1.1};
        std::mt19937_64 rng(16031605);
        std::size_t failures=0, tests=0;
        volatile double checksum=0.;
        std::cout<<std::setprecision(12);
        for(const auto n:{std::size_t(16),std::size_t(64),std::size_t(128)}) {
            auto a=sample(rng,n),b=sample(rng,n);
            const auto build_start=Clock::now();
            SpatialPersistentTree ta(a),tb(b);
            const auto build_end=Clock::now();
            const double root=local_spatial_aggregate_bound_guarded(ta.node(0).group,tb.node(0).group).value;
            for(const auto depth:depths) {
                const auto t0=Clock::now();
                const auto full=local_spatial_persistent_bound(ta,tb,depth);
                const auto t1=Clock::now();
                const auto recursive=local_spatial_adaptive_bound(a,b,depth);
                if(!close(full.value,recursive.value) || full.statistics.node_visits!=recursive.statistics.node_visits) ++failures;
                for(const double fraction:fractions) {
                    const double budget=root*fraction;
                    const auto start=Clock::now();
                    const auto stopped=local_spatial_persistent_bound(ta,tb,depth,budget,true);
                    const auto finish=Clock::now();
                    ++tests;
                    // A positive decision must always be supported by a valid bound.
                    // A negative decision may be conservative due to early termination.
                    if(stopped.within_budget && (!std::isfinite(stopped.value) || stopped.value>budget || full.value>budget)) ++failures;
                    if(stopped.statistics.node_visits>full.statistics.node_visits) ++failures;
                    if(fraction>=1. && !stopped.within_budget) ++failures;
                    if(fraction>=1. && stopped.statistics.node_visits!=1) ++failures;
                    checksum=checksum+stopped.value*1e-100;
                    std::cout<<"size="<<n<<" depth="<<depth<<" budget_fraction="<<fraction
                        <<" root_bound="<<root<<" full_bound="<<full.value
                        <<" full_nodes="<<full.statistics.node_visits
                        <<" budget_nodes="<<stopped.statistics.node_visits
                        <<" full_seconds="<<seconds(t0,t1)
                        <<" budget_seconds="<<seconds(start,finish)
                        <<" within_budget="<<stopped.within_budget
                        <<" failures="<<failures<<'\n';
                }
            }
            std::cout<<"size="<<n<<" tree_nodes_a="<<ta.node_count()
                     <<" tree_nodes_b="<<tb.node_count()
                     <<" build_seconds="<<seconds(build_start,build_end)<<'\n';
        }
        std::cout<<"cases="<<tests<<" failures="<<failures<<" checksum="<<double(checksum)
                 <<" formally_certified=0 "<<(failures?"FAIL":"OK")<<'\n';
        return failures?1:0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
