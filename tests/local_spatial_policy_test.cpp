#include "chem/kumo/local_spatial_decision.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace chem::kumo;
static void check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
int main() {
    try {
        const std::vector<GaussianPairEnvelope> a{{1.0,{0.,0.,0.},1.0}};
        const std::vector<GaussianPairEnvelope> b{{1.0,{1.,0.,0.},1.0}};
        SpatialPersistentTree ta(a),tb(b),empty(std::vector<GaussianPairEnvelope>{});
        const auto full=local_spatial_persistent_bound(ta,tb,8);
        check(std::isfinite(full.value)&&full.value>0.,"invalid reference");
        const double budget=full.value*2.;
        const auto strict=local_spatial_persistent_decision(ta,tb,8,budget,BoundPolicy::Rigorous);
        const auto experimental=local_spatial_persistent_decision(ta,tb,8,budget,BoundPolicy::Experimental);
        const auto legacy=local_spatial_persistent_decision(ta,tb,8,budget);
        check(strict.decision==SpatialDecision::KeepOrRefine,"rigorous accepted engineering bound");
        check(!strict.formally_certified,"incorrect certification");
        check(experimental.decision==SpatialDecision::CertifiedSkip,"experimental regression");
        check(legacy.decision==experimental.decision,"legacy regression");
        check(experimental.upper_bound==legacy.upper_bound,"legacy bound mismatch");
        for(double bad:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::infinity()}) {
            const auto s=local_spatial_persistent_decision(ta,tb,8,bad,BoundPolicy::Rigorous);
            check(s.decision==SpatialDecision::KeepOrRefine,"bad budget accepted");
        }
        const auto z=local_spatial_persistent_decision(empty,tb,8,0.,BoundPolicy::Rigorous);
        check(z.decision==SpatialDecision::CertifiedSkip && z.formally_certified && z.upper_bound==0.,
              "empty product should be exactly zero");
        std::cout<<"rigorous=keep experimental=skip legacy=compatible empty=certified status=OK\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
