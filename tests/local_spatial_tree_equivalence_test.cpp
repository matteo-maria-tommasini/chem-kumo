#include "chem/kumo/local_spatial_persistent_tree.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
using namespace chem::kumo;

static double reference(const std::vector<GaussianPairEnvelope>& a,
                        const std::vector<GaussianPairEnvelope>& b) {
    if(a.size()==1 && b.size()==1) return local_gaussian_envelope_overlap_bound(a[0],b[0]);
    long double s=0;
    for(const auto& x:a) for(const auto& y:b) s+=local_gaussian_envelope_overlap_bound(x,y);
    if(s>=std::numeric_limits<double>::max()) return INFINITY;
    return std::nextafter(static_cast<double>(s),INFINITY);
}
static bool check_partition(const SpatialPersistentTree& tree,std::size_t id,
                            std::vector<GaussianPairEnvelope> expected) {
    const auto& n=tree.node(id);
    if(n.end-n.begin!=expected.size()) return false;
    if(n.leaf) return expected.size()==1;
    auto [left,right]=spatial_adaptive_detail::split(expected);
    if(tree.node(n.left).begin!=n.begin || tree.node(n.left).end!=n.begin+left.size() ||
       tree.node(n.right).begin!=n.begin+left.size() || tree.node(n.right).end!=n.end) return false;
    auto same_members=[&](std::size_t begin,const std::vector<GaussianPairEnvelope>& items) {
        std::vector<bool> used(items.size(),false);
        for(std::size_t k=0;k<items.size();++k) {
            const auto& x=tree.leaf(begin+k);
            bool found=false;
            for(std::size_t j=0;j<items.size();++j)
                if(!used[j] && x.exponent==items[j].exponent && x.amplitude==items[j].amplitude && x.center==items[j].center) {
                    used[j]=true;found=true;break;
                }
            if(!found) return false;
        }
        return true;
    };
    if(!same_members(n.begin,left)||!same_members(n.begin+left.size(),right)) return false;
    return check_partition(tree,n.left,std::move(left)) && check_partition(tree,n.right,std::move(right));
}
int main() {
    std::mt19937_64 rng(160316);
    std::uniform_real_distribution<double> coordinate(-12,12),exponent(.2,8),amplitude(.01,2);
    std::size_t failures=0,partition_failures=0,numeric_failures=0,visits_mismatch=0;
    long double ratio_sum=0;
    double worst_scaled=0;
    for(int trial=0;trial<500;++trial) {
        std::vector<GaussianPairEnvelope> a,b;
        const int na=1+trial%13,nb=1+(trial*7)%17;
        for(int group=0;group<2;++group) {
            auto& out=group==0?a:b;
            for(int i=0;i<(group==0?na:nb);++i) {
                GaussianPairEnvelope e;
                e.exponent=exponent(rng);e.amplitude=amplitude(rng);
                for(double& c:e.center)c=coordinate(rng)+(group==0?0.0:5.0);
                out.push_back(e);
            }
        }
        SpatialPersistentTree ta(a),tb(b);
        const bool partitions=check_partition(ta,0,a)&&check_partition(tb,0,b);
        const auto old=local_spatial_adaptive_bound(a,b,8);
        const auto now=local_spatial_persistent_bound(ta,tb,8);
        const double ref=reference(a,b);
        const double scaled=std::abs(old.value-now.value)/std::max({1e-300,std::abs(old.value),std::abs(now.value)});
        const bool visits=old.statistics.node_visits==now.statistics.node_visits &&
            old.statistics.splits==now.statistics.splits && old.statistics.leaf_pairs==now.statistics.leaf_pairs;
        const bool valid=now.value>=ref && old.value>=ref && scaled<=1e-13;
        if(!partitions) ++partition_failures;
        if(!visits) ++visits_mismatch;
        if(!valid) ++numeric_failures;
        worst_scaled=std::max(worst_scaled,scaled);
        if(ref>0 && std::isfinite(ref)) ratio_sum+=now.value/ref;
        if((!partitions||!visits||!valid) && failures++==0) {
            std::cerr<<std::setprecision(17)<<"FAIL trial="<<trial<<" na="<<na<<" nb="<<nb
                     <<" partitions="<<partitions<<" visits_equal="<<visits
                     <<" reference="<<ref<<" recursive="<<old.value<<" persistent="<<now.value
                     <<" scaled_difference="<<scaled<<" old_nodes="<<old.statistics.node_visits
                     <<" new_nodes="<<now.statistics.node_visits<<'\n';
        }
    }
    std::cout<<std::setprecision(12)<<"cases=500 failures="<<failures
             <<" partition_failures="<<partition_failures<<" numeric_failures="<<numeric_failures
             <<" visits_mismatch="<<visits_mismatch<<" worst_scaled_difference="<<worst_scaled
             <<" mean_persistent_to_detailed="<<double(ratio_sum/500)
             <<" formally_certified=0 "<<(failures?"FAIL":"OK")<<'\n';
    return failures?1:0;
}
