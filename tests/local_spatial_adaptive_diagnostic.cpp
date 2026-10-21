#include "chem/kumo/local_spatial_adaptive_bound.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;

static double detailed(const std::vector<GaussianPairEnvelope>& a,
                       const std::vector<GaussianPairEnvelope>& b) {
    if(a.size()==1 && b.size()==1) return local_gaussian_envelope_overlap_bound(a.front(),b.front());
    long double sum=0;
    for(const auto& x:a) for(const auto& y:b)
        sum+=local_gaussian_envelope_overlap_bound(x,y);
    if(sum>=std::numeric_limits<double>::max()) return INFINITY;
    return std::nextafter(static_cast<double>(sum),INFINITY);
}
int main() {
    std::mt19937_64 rng(160316);
    std::uniform_real_distribution<double> coordinate(-12,12), exponent(0.2,8), amplitude(0.01,2);
    std::size_t cases=0,failures=0,improvements=0;
    double max_relative_shortfall=0;
    long double mean_root_ratio=0,mean_adaptive_ratio=0;
    std::size_t total_nodes=0;
    for(int trial=0;trial<500;++trial) {
        std::vector<GaussianPairEnvelope> a,b;
        const int na=1+(trial%13),nb=1+((trial*7)%17);
        for(int group=0;group<2;++group) {
            auto& out=group==0?a:b;
            for(int i=0;i<(group==0?na:nb);++i) {
                GaussianPairEnvelope e;
                e.exponent=exponent(rng);e.amplitude=amplitude(rng);
                for(double& c:e.center)c=coordinate(rng)+(group==0?0.0:5.0);
                out.push_back(e);
            }
        }
        const double ref=detailed(a,b);
        const double root=local_spatial_aggregate_bound_guarded(
            make_spatial_envelope_group(a),make_spatial_envelope_group(b)).value;
        const auto adaptive=local_spatial_adaptive_bound(a,b,8);
        if(adaptive.value>root) throw std::runtime_error("adaptive bound exceeds root");
        if(adaptive.value<ref) {
            ++failures;
            max_relative_shortfall=std::max(max_relative_shortfall,(ref-adaptive.value)/ref);
            std::cerr<<std::setprecision(17)<<"FAIL trial="<<trial<<" na="<<na<<" nb="<<nb
                     <<" reference="<<ref<<" root="<<root<<" adaptive="<<adaptive.value<<'\n';
        }
        if(adaptive.value<root) ++improvements;
        if(ref>0 && std::isfinite(ref)) {
            mean_root_ratio+=root/ref;
            mean_adaptive_ratio+=adaptive.value/ref;
        }
        total_nodes+=adaptive.statistics.node_visits;
        ++cases;
    }
    std::cout<<std::setprecision(12)<<"cases="<<cases<<" failures="<<failures
             <<" improvements="<<improvements<<" mean_root_ratio="<<double(mean_root_ratio/cases)
             <<" mean_adaptive_ratio="<<double(mean_adaptive_ratio/cases)
             <<" mean_node_visits="<<double(total_nodes)/cases
             <<" max_relative_shortfall="<<max_relative_shortfall
             <<" formally_certified=0"<<(failures?" FAIL":" OK")<<'\n';
    return failures?1:0;
}
