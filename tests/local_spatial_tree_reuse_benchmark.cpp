#include "chem/kumo/local_spatial_persistent_tree.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;
using Clock = std::chrono::steady_clock;
static double elapsed(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b-a).count();
}
static std::vector<GaussianPairEnvelope> make_group(std::mt19937_64& rng, std::size_t n) {
    std::uniform_real_distribution<double> position(-12.,12.), exponent(.2,8.), amplitude(.01,2.);
    std::vector<GaussianPairEnvelope> result;
    result.reserve(n);
    for(std::size_t i=0;i<n;++i)
        result.push_back({exponent(rng),{position(rng),position(rng),position(rng)},amplitude(rng)});
    return result;
}
static bool close(double a,double b) {
    if(a==b) return true;
    if(!std::isfinite(a)||!std::isfinite(b)) return false;
    return std::abs(a-b)<=1e-12*std::max({1.,std::abs(a),std::abs(b)});
}
int main(int argc,char** argv) {
    try {
        if(argc>2) throw std::invalid_argument("usage: benchmark [--quick]");
        const bool quick=argc==2 && std::string(argv[1])=="--quick";
        if(argc==2 && !quick) throw std::invalid_argument("unknown option");
        const std::vector<std::size_t> sizes=quick ? std::vector<std::size_t>{16,32} : std::vector<std::size_t>{16,32,64,128};
        const std::vector<std::size_t> repeats=quick ? std::vector<std::size_t>{1,10} : std::vector<std::size_t>{1,10,100,1000};
        std::mt19937_64 rng(16031604);
        volatile double sink=0.;
        std::size_t failures=0;
        std::cout<<std::setprecision(12);
        for(auto n:sizes) {
            auto a=make_group(rng,n), b=make_group(rng,n);
            const auto t0=Clock::now();
            SpatialPersistentTree ta(a),tb(b);
            const auto t1=Clock::now();
            const double build=elapsed(t0,t1);
            // Verify identical numerical work on exactly the same groups.
            const auto reference=local_spatial_adaptive_bound(a,b,8);
            const auto persistent=local_spatial_persistent_bound(ta,tb,8);
            if(!close(reference.value,persistent.value) ||
               reference.statistics.node_visits!=persistent.statistics.node_visits) ++failures;
            // A budget at the root bound is deliberately chosen so the early-exit path is exercised.
            const double budget=local_spatial_aggregate_bound_guarded(
                ta.node(0).group,tb.node(0).group).value;
            const auto stopped=local_spatial_persistent_bound(ta,tb,8,budget,true);
            if(!stopped.within_budget || !close(stopped.value,budget)) ++failures;
            for(auto count:repeats) {
                auto start=Clock::now();
                for(std::size_t i=0;i<count;++i) {
                    auto r=local_spatial_adaptive_bound(a,b,8);
                    sink=sink+r.value*1e-100;
                }
                auto mid=Clock::now();
                for(std::size_t i=0;i<count;++i) {
                    auto r=local_spatial_persistent_bound(ta,tb,8);
                    sink=sink+r.value*1e-100;
                }
                auto end=Clock::now();
                for(std::size_t i=0;i<count;++i) {
                    auto r=local_spatial_persistent_bound(ta,tb,8,budget,true);
                    sink=sink+r.value*1e-100;
                }
                auto finish=Clock::now();
                const double recursive=elapsed(start,mid), persistent_time=elapsed(mid,end);
                const double budget_time=elapsed(end,finish);
                const double recursive_per=recursive/double(count), persistent_per=persistent_time/double(count);
                const double saving=recursive_per-persistent_per;
                const double break_even=saving>0 ? std::ceil(build/saving) : std::numeric_limits<double>::infinity();
                std::cout<<"size="<<n<<" queries="<<count
                    <<" build_seconds="<<build<<" recursive_seconds="<<recursive
                    <<" persistent_seconds="<<persistent_time<<" budget_seconds="<<budget_time
                    <<" recursive_nodes="<<reference.statistics.node_visits
                    <<" persistent_nodes="<<persistent.statistics.node_visits
                    <<" budget_nodes="<<stopped.statistics.node_visits
                    <<" speedup_query="<<(persistent_time>0 ? recursive/persistent_time : 0.)
                    <<" break_even_queries="<<break_even
                    <<" failures="<<failures<<'\n';
            }
        }
        std::cout<<"failures="<<failures<<" checksum="<<double(sink)
                 <<" formally_certified=0 "<<(failures?"FAIL":"OK")<<'\n';
        return failures?1:0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
