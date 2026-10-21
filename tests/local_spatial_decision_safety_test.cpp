#include "chem/kumo/local_spatial_decision.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
static void require(bool ok,const char* msg) {if(!ok) throw std::runtime_error(msg);}
int main() {
    try {
        const std::vector<GaussianPairEnvelope> a{{1.0,{0.0,0.0,0.0},1.0},
                                                   {2.0,{1.0,0.0,0.0},0.5}};
        const std::vector<GaussianPairEnvelope> b{{1.5,{0.5,0.0,0.0},1.0},
                                                   {0.8,{2.0,0.0,0.0},0.25}};
        SpatialPersistentTree ta(a),tb(b),empty({});
        const auto full=local_spatial_persistent_bound(ta,tb,8);
        require(std::isfinite(full.value)&&full.value>0,"invalid reference");
        std::size_t checks=0, skips=0;
        for(double budget: {0.0,full.value*0.01,full.value*0.5,full.value,
                            full.value*2.0,-1.0,
                            std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::infinity()}) {
            const auto r=local_spatial_persistent_decision(ta,tb,8,budget);
            require(!r.formally_certified,"unexpected formal certification");
            if(r.decision==SpatialDecision::CertifiedSkip) {
                require(std::isfinite(budget)&&budget>0&&r.upper_bound<=budget,
                        "unsafe positive decision");
                require(full.value<=budget,"positive decision contradicts full bound");
                ++skips;
            }
            if(!std::isfinite(budget)||budget<=0)
                require(r.decision==SpatialDecision::KeepOrRefine,"invalid budget accepted");
            ++checks;
        }
        for(double budget: {0.0,1.0}) {
            const auto r=local_spatial_persistent_decision(empty,tb,8,budget);
            require(r.decision==SpatialDecision::CertifiedSkip && r.upper_bound==0,
                    "empty group should be zero");
            ++checks;
        }
        const auto strict=local_spatial_persistent_decision(ta,tb,8,0.0);
        require(strict.decision==SpatialDecision::KeepOrRefine,"zero budget must keep");
        const auto easy=local_spatial_persistent_decision(ta,tb,8,full.value*2);
        require(easy.decision==SpatialDecision::CertifiedSkip,"easy budget must skip");
        std::cout<<"checks="<<checks<<" skips="<<skips
                 <<" zero_budget=keep invalid_budgets=keep empty=skip"
                 <<" formally_certified=0 OK\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
