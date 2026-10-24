#include "chem/kumo/fchk.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    using namespace chem::kumo;
    if (argc < 2) { std::cerr << "usage: diagnostic file.fchk[.gz] ...\n"; return 2; }
    try {
        constexpr std::array<double,4> relative_thresholds{1e-2,1e-4,1e-6,1e-8};
        for (int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells=local_shell_components(basis);
            const auto norms=local_ao_product_norms(basis);
            const double reference_scale=norms.maxCoeff()*norms.maxCoeff();
            std::vector<Eigen::Index> offsets(shells.size()+1,0);
            for(std::size_t a=0;a<shells.size();++a)
                offsets[a+1]=offsets[a]+static_cast<Eigen::Index>(shells[a].size());
            // Centers from diagonal shell-pair envelopes. This avoids assumptions
            // about the internal primitive representation.
            std::vector<std::array<double,3>> centers;
            centers.reserve(shells.size());
            for(std::size_t a=0;a<shells.size();++a) {
                const auto e=prepare_local_shell_pair_envelope(shells[a],shells[a],0.7);
                if(e.at(0,0).empty()) throw std::runtime_error("empty diagonal shell envelope");
                centers.push_back(e.at(0,0).front().center);
            }
            struct Pick {std::size_t a,b; double distance;};
            std::array<std::vector<Pick>,4> picks;
            constexpr std::size_t max_per_bin=8;
            // Deterministic, distributed shell sampling: one representative
            // per distance stratum from each part of the molecule.
            for(std::size_t a=0;a<shells.size();++a)
                for(std::size_t b=a+1;b<shells.size();++b) {
                    double d2=0;
                    for(int t=0;t<3;++t) d2+=std::pow(centers[a][t]-centers[b][t],2);
                    const double d=std::sqrt(d2);
                    const std::size_t bin=d<2?0:d<5?1:d<10?2:3;
                    if(picks[bin].size()<max_per_bin) picks[bin].push_back({a,b,d});
                    else {
                        // Reservoir-like deterministic spread across shell indices.
                        const std::size_t slot=(a*1315423911ULL+b*2654435761ULL)%max_per_bin;
                        if ((a+b)%17==0) picks[bin][slot]={a,b,d};
                    }
                }
            std::cout<<std::scientific<<std::setprecision(10);
            for(std::size_t bin=0;bin<4;++bin) {
                std::size_t elements=0,tighter=0,violations=0;
                std::array<std::size_t,4> cs_skips{},hybrid_skips{};
                long double logratio=0,max_scaled_violation=0;
                double min_distance=std::numeric_limits<double>::infinity(),max_distance=0;
                for(const auto& pick:picks[bin]) {
                    const auto a=pick.a,b=pick.b;
                    min_distance=std::min(min_distance,pick.distance);
                    max_distance=std::max(max_distance,pick.distance);
                    // Compare (aa|bb): product centers are separated by
                    // the shell-center distance. This isolates spatial decay.
                    const auto aa=prepare_local_shell_pair_envelope(shells[a],shells[a],0.7);
                    const auto bb=prepare_local_shell_pair_envelope(shells[b],shells[b],0.7);
                    const auto q=local_shell_quartet({&shells[a],&shells[a],&shells[b],&shells[b]});
                    for(std::size_t i=0;i<aa.na;++i)
                      for(std::size_t j=0;j<aa.nb;++j)
                        for(std::size_t k=0;k<bb.na;++k)
                          for(std::size_t l=0;l<bb.nb;++l) {
                            const auto mu=offsets[a]+static_cast<Eigen::Index>(i);
                            const auto nu=offsets[a]+static_cast<Eigen::Index>(j);
                            const auto la=offsets[b]+static_cast<Eigen::Index>(k);
                            const auto si=offsets[b]+static_cast<Eigen::Index>(l);
                            const double cs=norms(mu,nu)*norms(la,si);
                            const double sp=local_shell_pair_envelope_bound(aa,i,j,bb,k,l);
                            if(!std::isfinite(cs)||!std::isfinite(sp)||cs<0||sp<0)
                                throw std::runtime_error("invalid bound");
                            const double hy=std::min(cs,sp);
                            const double exact=std::abs(q.at(i,j,k,l));
                            const long double scale=std::max({1.L,(long double)exact,(long double)hy});
                            const long double excess=std::max(0.L,(long double)exact-hy);
                            max_scaled_violation=std::max(max_scaled_violation,excess/scale);
                            if(excess>2e-11L*scale) ++violations;
                            tighter+=(sp<cs);
                            if(cs>0 && hy>0) logratio+=std::log((long double)hy/cs);
                            // Threshold relative to the global maximum CS scale.
                            // This is a diagnostic, not a cumulative error budget.
                            for(std::size_t t=0;t<relative_thresholds.size();++t) {
                                const double threshold=relative_thresholds[t]*reference_scale;
                                cs_skips[t]+=(cs<threshold);
                                hybrid_skips[t]+=(hy<threshold);
                            }
                            ++elements;
                          }
                }
                std::cout<<"system="<<argv[arg]<<" bin="<<bin
                         <<" range="<<(bin==0?"0-2":bin==1?"2-5":bin==2?"5-10":"10+")
                         <<" shell_pairs="<<picks[bin].size()
                         <<" distance_min="<<(picks[bin].empty()?0:min_distance)
                         <<" distance_max="<<max_distance
                         <<" ao_elements="<<elements
                         <<" spatial_tighter="<<tighter
                         <<" tighter_fraction="<<(elements?double(tighter)/elements:0)
                         <<" geometric_mean_hybrid_to_cs="<<(elements?std::exp(logratio/elements):1.L)
                         <<" max_scaled_violation="<<max_scaled_violation
                         <<" violations="<<violations;
                for(std::size_t t=0;t<relative_thresholds.size();++t)
                    std::cout<<" relative_threshold="<<relative_thresholds[t]
                             <<" cs_skips="<<cs_skips[t]
                             <<" hybrid_skips="<<hybrid_skips[t];
                std::cout<<" status="<<(violations?"FAIL":"OK")<<'\n';
                if(violations) return 1;
            }
        }
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
