#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/local_spatial_aggregate_bound.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using chem::kumo::GaussianPairEnvelope;
struct Group {
    std::array<long double,3> center{};
    long double radius=0, min_exponent=std::numeric_limits<long double>::infinity();
    long double amplitude_sum=0;
    std::vector<GaussianPairEnvelope> leaves;
};
Group make_group(const std::vector<GaussianPairEnvelope>& leaves) {
    Group g;
    g.leaves=leaves;
    if (leaves.empty()) return g;
    for (const auto& e: leaves) {
        if (!(e.exponent>0) || !(e.amplitude>=0) || !std::isfinite(e.exponent) || !std::isfinite(e.amplitude))
            throw std::runtime_error("nonfinite envelope; cannot construct group");
        for (int d=0;d<3;++d) g.center[d]+=static_cast<long double>(e.center[d])/leaves.size();
        g.min_exponent=std::min(g.min_exponent,static_cast<long double>(e.exponent));
        g.amplitude_sum+=static_cast<long double>(e.amplitude);
    }
    for (const auto& e: leaves) {
        long double r2=0;
        for(int d=0;d<3;++d) {
            const long double x=static_cast<long double>(e.center[d])-g.center[d];
            r2+=x*x;
        }
        g.radius=std::max(g.radius,std::sqrt(r2));
    }
    return g;
}
long double group_bound(const Group& a,const Group& b) {
    if(a.leaves.empty()||b.leaves.empty()) return 0;
    long double d2=0;
    for(int k=0;k<3;++k) {
        const long double v=a.center[k]-b.center[k]; d2+=v*v;
    }
    const long double separation=std::max(0.0L,std::sqrt(d2)-a.radius-b.radius);
    const long double p=a.min_exponent+b.min_exponent;
    const long double eta=a.min_exponent*b.min_exponent/p;
    const long double pi=std::acos(-1.0L);
    return a.amplitude_sum*b.amplitude_sum*std::pow(pi/p,1.5L)*std::exp(-eta*separation*separation);
}
long double detailed_bound(const Group& a,const Group& b) {
    long double total=0;
    for(const auto& x:a.leaves) for(const auto& y:b.leaves)
        total+=chem::kumo::local_gaussian_envelope_overlap_bound(x,y);
    return total;
}
long double detailed_extended(const Group& a,const Group& b) {
    const long double pi=std::acos(-1.0L);
    long double total=0;
    for(const auto& x:a.leaves) for(const auto& y:b.leaves) {
        const long double alpha=x.exponent, beta=y.exponent;
        const long double p=alpha+beta, eta=alpha*beta/p;
        long double r2=0;
        for(int d=0;d<3;++d) {
            const long double delta=static_cast<long double>(x.center[d])-y.center[d];
            r2+=delta*delta;
        }
        total+=static_cast<long double>(x.amplitude)*y.amplitude*
            std::pow(pi/p,1.5L)*std::exp(-eta*r2);
    }
    return total;
}
}
int main(int argc,char** argv) {
    using namespace chem::kumo;
    if(argc<2) {std::cerr<<"usage: diagnostic file.fchk.gz [file2...]\n";return 2;}
    try {
        for(int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells=local_shell_components(basis);
            const std::size_t ns=std::min<std::size_t>(6,shells.size());
            std::vector<Group> groups;
            // Use a bounded sample, deliberately independent of the number of HBC shell pairs.
            for(std::size_t a=0;a<ns;++a) for(std::size_t b=a;b<ns;++b) {
                const auto e=prepare_local_shell_pair_envelope(shells[a],shells[b],0.7);
                std::vector<GaussianPairEnvelope> leaves;
                for(const auto& component:e.components)
                    leaves.insert(leaves.end(),component.begin(),component.end());
                groups.push_back(make_group(leaves));
            }
            std::size_t checked=0,tighter=0, extended_violations=0, api_violations=0;
            std::size_t guarded_violations=0, guarded_nonfinite=0, singleton_pairs=0;
            long double max_guarded_shortfall=0, max_guarded_relative_shortfall=0;
            long double guarded_ratio_sum=0, guarded_max_ratio=0;
            long double max_violation=0,ratio_sum=0, max_ratio=0;
            long double max_relative_extended=0, max_relative_api=0;
            struct Finding { long double relative, absolute, aggregate, extended, api; std::size_t i,j,na,nb; };
            std::vector<Finding> findings;
            auto start=std::chrono::steady_clock::now();
            for(std::size_t i=0;i<groups.size();++i) for(std::size_t j=i;j<groups.size();++j) {
                const long double aggregate=group_bound(groups[i],groups[j]);
                const long double detailed=detailed_bound(groups[i],groups[j]);
                const long double extended=detailed_extended(groups[i],groups[j]);
                const auto guarded=local_spatial_aggregate_bound_guarded(
                    make_spatial_envelope_group(groups[i].leaves),
                    make_spatial_envelope_group(groups[j].leaves));
                if(groups[i].leaves.size()==1 && groups[j].leaves.size()==1) ++singleton_pairs;
                const long double gv=static_cast<long double>(guarded.value);
                if(!std::isfinite(gv)) ++guarded_nonfinite;
                if(gv<detailed) {
                    ++guarded_violations;
                    max_guarded_shortfall=std::max(max_guarded_shortfall,detailed-gv);
                    if(detailed>0) max_guarded_relative_shortfall=std::max(
                        max_guarded_relative_shortfall,(detailed-gv)/detailed);
                    std::cerr<<std::scientific<<std::setprecision(18)
                        <<"guarded_failure system="<<argv[arg]<<" block_i="<<i
                        <<" block_j="<<j<<" guarded="<<gv<<" detailed="<<detailed
                        <<" shortfall="<<(detailed-gv)<<"\n";
                }
                if(detailed>0) {
                    const long double gr=gv/detailed;
                    guarded_ratio_sum+=gr;
                    guarded_max_ratio=std::max(guarded_max_ratio,gr);
                }
                if(!std::isfinite(aggregate)||!std::isfinite(detailed)) continue;
                max_violation=std::max(max_violation,detailed-aggregate);
                const long double scale=std::max(std::abs(extended),std::abs(aggregate));
                const long double rel_extended=scale>0 ? std::max(0.0L,(extended-aggregate)/scale):0;
                const long double rel_api=std::max(std::abs(detailed),std::abs(aggregate))>0 ?
                    std::max(0.0L,(detailed-aggregate)/std::max(std::abs(detailed),std::abs(aggregate))):0;
                max_relative_extended=std::max(max_relative_extended,rel_extended);
                max_relative_api=std::max(max_relative_api,rel_api);
                if(extended>aggregate) ++extended_violations;
                if(detailed>aggregate) {
                    ++api_violations;
                    findings.push_back({rel_api,detailed-aggregate,aggregate,extended,detailed,i,j,
                        groups[i].leaves.size(),groups[j].leaves.size()});
                }
                if(aggregate+1e-12L*std::max(1.0L,detailed)<detailed)
                    throw std::runtime_error("aggregate bound is smaller than detailed bound");
                if(aggregate<detailed) ++tighter;
                if(detailed>0) {
                    const long double ratio=aggregate/detailed;
                    ratio_sum+=ratio;
                    max_ratio=std::max(max_ratio,ratio);
                }
                ++checked;
            }
            std::sort(findings.begin(),findings.end(),[](const Finding& a,const Finding& b) {
                return a.relative>b.relative;
            });
            for(std::size_t z=0;z<std::min<std::size_t>(10,findings.size());++z) {
                const auto& f=findings[z];
                std::cout<<std::scientific<<std::setprecision(18)
                    <<"audit_rank="<<(z+1)<<" block_i="<<f.i<<" block_j="<<f.j
                    <<" leaves_i="<<f.na<<" leaves_j="<<f.nb
                    <<" aggregate="<<f.aggregate<<" detailed_extended="<<f.extended
                    <<" detailed_api="<<f.api<<" absolute_violation="<<f.absolute
                    <<" relative_violation="<<f.relative<<"\n";
            }
            const auto stop=std::chrono::steady_clock::now();
            std::cout<<std::scientific<<std::setprecision(12)
                <<"system="<<argv[arg]<<" sampled_shells="<<ns
                <<" sampled_shell_pairs="<<groups.size()<<" block_comparisons="<<checked
                <<" max_positive_violation="<<max_violation
                <<" api_violations="<<api_violations
                <<" guarded_violations="<<guarded_violations
                <<" singleton_pairs="<<singleton_pairs
                <<" guarded_nonfinite="<<guarded_nonfinite
                <<" max_guarded_shortfall="<<max_guarded_shortfall
                <<" max_guarded_relative_shortfall="<<max_guarded_relative_shortfall
                <<" mean_guarded_to_detailed="<<(checked?guarded_ratio_sum/checked:0)
                <<" max_guarded_to_detailed="<<guarded_max_ratio
                <<" extended_violations="<<extended_violations
                <<" max_relative_api="<<max_relative_api
                <<" max_relative_extended="<<max_relative_extended
                <<" aggregate_tighter_due_roundoff="<<tighter
                <<" mean_aggregate_to_detailed="<<(checked?ratio_sum/checked:0)
                <<" max_aggregate_to_detailed="<<max_ratio
                <<" comparison_seconds="<<std::chrono::duration<double>(stop-start).count()
                <<" status="<<(guarded_violations==0?"OK":"FAIL")<<"\n";
            if(guarded_violations) throw std::runtime_error("guarded aggregate underestimates detailed API bound");
        }
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
