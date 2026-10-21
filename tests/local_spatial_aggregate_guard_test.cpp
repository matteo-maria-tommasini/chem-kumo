#include "chem/kumo/local_spatial_aggregate_bound.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;

int main() {
    try {
        std::mt19937_64 rng(1603);
        std::uniform_real_distribution<double> coordinate(-12,12), exponent(0.05,20), amplitude(0.001,3);
        // Reproduce the original failure at trial 64 before random testing.
        {
            GaussianPairEnvelope x{2.1536602962077235901,
                {2.43325856375162263845,-3.70962147723840418578,-2.24150689752122511678},
                2.15454487149636664256};
            GaussianPairEnvelope y{16.2609898352078729999,
                {-1.07502369381155027384,3.03005110493031626362,11.3121256450182485764},
                2.04563843862058014622};
            const auto r=local_spatial_aggregate_bound_guarded(
                make_spatial_envelope_group({x}),make_spatial_envelope_group({y}));
            if(r.value!=local_gaussian_envelope_overlap_bound(x,y))
                throw std::runtime_error("singleton API mismatch");
        }
        std::size_t checked=0;
        long double max_ratio=0,max_guard_fraction=0;
        for(std::size_t trial=0;trial<4000;++trial) {
            const std::size_t na=1+(trial%32),nb=1+((trial*7)%32);
            std::vector<GaussianPairEnvelope> a,b;
            for(std::size_t i=0;i<na+nb;++i) {
                GaussianPairEnvelope e;
                e.exponent=exponent(rng);e.amplitude=amplitude(rng);
                for(int d=0;d<3;++d)e.center[d]=coordinate(rng);
                (i<na?a:b).push_back(e);
            }
            const auto ga=make_spatial_envelope_group(a);
            const auto gb=make_spatial_envelope_group(b);
            const auto bound=local_spatial_aggregate_bound_guarded(ga,gb);
            long double detailed=0;
            for(const auto& x:a)for(const auto& y:b)
                detailed+=local_gaussian_envelope_overlap_bound(x,y);
            if(!(static_cast<long double>(bound.value)>=detailed)) {
                const long double excess=detailed-static_cast<long double>(bound.value);
                const auto p=ga.minimum_exponent+gb.minimum_exponent;
                const auto eta=ga.minimum_exponent*gb.minimum_exponent/p;
                long double center_distance2=0;
                for(int axis=0;axis<3;++axis) {
                    const auto dx=ga.center[axis]-gb.center[axis];
                    center_distance2+=dx*dx;
                }
                const auto separation=std::max(0.0L,std::sqrt(center_distance2)-ga.radius-gb.radius);
                std::cerr << std::setprecision(21)
                    << " exponent_audit trial=" << trial
                    << " eta=" << eta
                    << " center_distance2=" << center_distance2
                    << " separation=" << separation
                    << " exponent_argument=" << -eta*separation*separation
                    << " log_aggregate=" << (bound.unguarded>0 ? std::log(bound.unguarded) : -std::numeric_limits<long double>::infinity())
                    << " log_detailed=" << (detailed>0 ? std::log(detailed) : -std::numeric_limits<long double>::infinity()) << "\n";
                std::cerr << std::setprecision(21)
                    << "FAIL trial=" << trial << " na=" << na << " nb=" << nb
                    << " aggregate=" << bound.value
                    << " raw=" << bound.unguarded
                    << " guard=" << bound.guard
                    << " detailed=" << detailed
                    << " absolute_excess=" << excess
                    << " relative_excess=" << (detailed>0 ? excess/detailed : 0)
                    << " alpha_a=" << ga.minimum_exponent
                    << " alpha_b=" << gb.minimum_exponent
                    << " amplitude_a=" << ga.amplitude_sum
                    << " amplitude_b=" << gb.amplitude_sum
                    << " radius_a=" << ga.radius
                    << " radius_b=" << gb.radius << '\n';
                for(std::size_t i=0;i<a.size();++i)
                    std::cerr << "A " << i << " " << a[i].exponent << " "
                        << a[i].amplitude << " " << a[i].center[0] << " "
                        << a[i].center[1] << " " << a[i].center[2] << '\n';
                for(std::size_t i=0;i<b.size();++i)
                    std::cerr << "B " << i << " " << b[i].exponent << " "
                        << b[i].amplitude << " " << b[i].center[0] << " "
                        << b[i].center[1] << " " << b[i].center[2] << '\n';
                throw std::runtime_error("aggregate underestimates detailed bound");
            }
            if(bound.formally_certified)throw std::runtime_error("incorrect certification flag");
            if(detailed>0)max_ratio=std::max(max_ratio,static_cast<long double>(bound.value)/detailed);
            if(bound.unguarded>0)max_guard_fraction=std::max(max_guard_fraction,bound.guard/bound.unguarded);
            ++checked;
        }
        std::cout<<std::setprecision(12)<<"random_cases="<<checked
                 <<" max_ratio="<<max_ratio<<" max_guard_fraction="<<max_guard_fraction
                 <<" formally_certified=0 OK\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
