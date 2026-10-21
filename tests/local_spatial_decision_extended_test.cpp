#include "chem/kumo/local_spatial_decision.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>
using namespace chem::kumo;
static std::vector<GaussianPairEnvelope> sample(std::mt19937_64& rng,std::size_t n) {
    std::uniform_real_distribution<double> xyz(-12.,12.), alpha(.2,8.), amp(.01,2.);
    std::vector<GaussianPairEnvelope> v;
    for(std::size_t i=0;i<n;++i) v.push_back({alpha(rng),{xyz(rng),xyz(rng),xyz(rng)},amp(rng)});
    return v;
}
static void require(bool condition,const char* message) {
    if(!condition) throw std::runtime_error(message);
}
int main() {
    try {
        constexpr std::array<std::size_t,6> depths{2,4,6,8,10,12};
        constexpr std::array<double,8> fractions{0.,.01,.1,.5,.9,.99,1.,1.1};
        std::mt19937_64 rng(16031605);
        std::size_t cases=0, skips=0, keeps=0, conservative_keeps=0, early_exits=0;
        for(const auto n:{std::size_t(16),std::size_t(64),std::size_t(128)}) {
            const SpatialPersistentTree a(sample(rng,n)),b(sample(rng,n));
            const double root=local_spatial_aggregate_bound_guarded(a.node(0).group,b.node(0).group).value;
            require(std::isfinite(root)&&root>0.,"invalid root bound");
            for(const auto depth:depths) {
                const auto full=local_spatial_persistent_bound(a,b,depth);
                require(std::isfinite(full.value)&&full.value>=0.,"invalid full bound");
                for(const double fraction:fractions) {
                    const double budget=root*fraction;
                    const auto result=local_spatial_persistent_decision(a,b,depth,budget);
                    ++cases;
                    require(!result.formally_certified,"unexpected formal certification");
                    if(result.decision==SpatialDecision::CertifiedSkip) {
                        ++skips;
                        require(budget>0. && std::isfinite(budget),"invalid budget accepted");
                        require(std::isfinite(result.upper_bound)&&result.upper_bound>=0.,"invalid accepted bound");
                        require(result.upper_bound<=budget,"accepted bound exceeds budget");
                        require(full.value<=budget,"false positive against full query");
                    } else {
                        ++keeps;
                        if(full.value<=budget) ++conservative_keeps;
                    }
                    if(result.node_visits<full.statistics.node_visits) ++early_exits;
                    require(result.node_visits<=full.statistics.node_visits,"too many decision visits");
                    if(fraction==0.) require(result.decision==SpatialDecision::KeepOrRefine,"zero budget accepted");
                    if(fraction>=1.) require(result.decision==SpatialDecision::CertifiedSkip && result.node_visits==1,"easy budget rejected");
                }
                for(const double invalid:{-1.,std::numeric_limits<double>::quiet_NaN(),
                                          std::numeric_limits<double>::infinity(),
                                          -std::numeric_limits<double>::infinity()}) {
                    const auto result=local_spatial_persistent_decision(a,b,depth,invalid);
                    require(result.decision==SpatialDecision::KeepOrRefine && result.node_visits==0,
                            "invalid budget must not query tree");
                }
            }
        }
        require(cases==144,"incorrect sweep size");
        const SpatialPersistentTree empty(std::vector<GaussianPairEnvelope>{});
        const SpatialPersistentTree nonempty(std::vector<GaussianPairEnvelope>{{1.,{0.,0.,0.},1.}});
        for(double budget:{0.,1.}) {
            const auto r=local_spatial_persistent_decision(empty,nonempty,8,budget);
            require(r.decision==SpatialDecision::CertifiedSkip && r.upper_bound==0.,"empty input not skipped");
        }
        std::cout<<"cases="<<cases<<" skips="<<skips<<" keeps="<<keeps
                 <<" conservative_keeps="<<conservative_keeps<<" early_exits="<<early_exits
                 <<" invalid_budgets=keep zero_budget=keep empty=skip"
                 <<" formally_certified=0 OK\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
