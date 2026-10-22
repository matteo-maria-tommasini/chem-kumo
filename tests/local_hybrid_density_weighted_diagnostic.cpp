#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"
#include <algorithm>
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
        for (int arg=1; arg<argc; ++arg) {
            const auto molecule=Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells=local_shell_components(basis);
            std::vector<Eigen::Index> offsets(shells.size()+1,0);
            for (std::size_t a=0;a<shells.size();++a)
                offsets[a+1]=offsets[a]+static_cast<Eigen::Index>(shells[a].size());
            const auto pt=molecule.density_matrix(Spin::Total);
            const auto ps=molecule.density_matrix(Spin::Difference);
            if (pt.rows()!=offsets.back() || ps.rows()!=offsets.back())
                throw std::runtime_error("density/basis dimension mismatch");
            const std::vector<Eigen::MatrixXd> densities{pt,ps};
            const auto norms=local_ao_product_norms(basis);
            // Spread the sample across the entire molecular basis; canonical shell pairs.
            const std::size_t ns=std::min<std::size_t>(6,shells.size());
            std::vector<std::size_t> chosen;
            for (std::size_t i=0;i<ns;++i) {
                const auto a=ns==1 ? 0 : i*(shells.size()-1)/(ns-1);
                if (chosen.empty() || chosen.back()!=a) chosen.push_back(a);
            }
            struct Pair { std::size_t a,b; LocalShellPairEnvelope env; };
            std::vector<Pair> pairs;
            for (std::size_t i=0;i<chosen.size();++i)
                for (std::size_t j=i;j<chosen.size();++j) {
                    const auto a=chosen[i], b=chosen[j];
                    pairs.push_back({a,b,prepare_local_shell_pair_envelope(shells[a],shells[b],0.7)});
                }
            // Bound the diagnostic cost for larger bases.
            constexpr std::size_t max_pairs=12;
            if (pairs.size()>max_pairs) pairs.resize(max_pairs);
            constexpr double thresholds[]{1e-6,1e-8,1e-10};
            std::size_t quartets=0, ao_elements=0, improved=0, failures=0;
            std::size_t cs_skips[3]{}, hybrid_skips[3]{};
            long double max_scaled_violation=0, max_cumulative_ratio=0;
            // For each destination AO product, accumulate absolute omitted-source
            // contributions. This is a sampled cumulative budget, NOT a global W bound.
            for (std::size_t dest=0;dest<pairs.size();++dest) {
                const auto& x=pairs[dest];
                for (std::size_t i=0;i<x.env.na;++i)
                    for (std::size_t j=0;j<x.env.nb;++j) {
                        const auto mu=offsets[x.a]+static_cast<Eigen::Index>(i);
                        const auto nu=offsets[x.b]+static_cast<Eigen::Index>(j);
                        long double cumulative_cs[2]{}, cumulative_hybrid[2]{}, cumulative_actual[2]{};
                        for (std::size_t src=0;src<pairs.size();++src) {
                            const auto& y=pairs[src];
                            const auto q=local_shell_quartet({&shells[x.a],&shells[x.b],&shells[y.a],&shells[y.b]});
                            ++quartets;
                            for (std::size_t k=0;k<y.env.na;++k)
                                for (std::size_t l=0;l<y.env.nb;++l) {
                                    const auto lambda=offsets[y.a]+static_cast<Eigen::Index>(k);
                                    const auto sigma=offsets[y.b]+static_cast<Eigen::Index>(l);
                                    const double cs=norms(mu,nu)*norms(lambda,sigma);
                                    const double spatial=local_shell_pair_envelope_bound(x.env,i,j,y.env,k,l);
                                    if (!std::isfinite(cs) || !std::isfinite(spatial) || cs<0 || spatial<0)
                                        throw std::runtime_error("invalid bound");
                                    const double hybrid=std::min(cs,spatial);
                                    const double actual=std::abs(q.at(i,j,k,l));
                                    const long double excess=std::max(0.0L,static_cast<long double>(actual)-hybrid);
                                    const long double scale=std::max({1.0L,static_cast<long double>(actual),static_cast<long double>(hybrid)});
                                    max_scaled_violation=std::max(max_scaled_violation,excess/scale);
                                    if (excess>2e-11L*scale) ++failures;
                                    improved+=(spatial<cs);
                                    for (std::size_t t=0;t<3;++t) {
                                        cs_skips[t]+=(cs<thresholds[t]);
                                        hybrid_skips[t]+=(hybrid<thresholds[t]);
                                    }
                                    // Canonical source shell pair: off-diagonal shells
                                    // have a reciprocal partner. For same shell, AO
                                    // off-diagonal components occur explicitly in loops.
                                    const long double m=y.a==y.b ? 1.0L : 2.0L;
                                    for (std::size_t ch=0;ch<2;++ch) {
                                        const long double w=m*std::abs(densities[ch](lambda,sigma));
                                        cumulative_cs[ch]+=w*cs;
                                        cumulative_hybrid[ch]+=w*hybrid;
                                        cumulative_actual[ch]+=w*actual;
                                    }
                                    ++ao_elements;
                                }
                        }
                        for (std::size_t ch=0;ch<2;++ch) {
                            if (cumulative_actual[ch]>cumulative_hybrid[ch]+2e-11L*std::max(1.0L,cumulative_actual[ch])) ++failures;
                            if (cumulative_cs[ch]>0)
                                max_cumulative_ratio=std::max(max_cumulative_ratio,cumulative_hybrid[ch]/cumulative_cs[ch]);
                        }
                    }
            }
            std::cout<<std::scientific<<std::setprecision(12)
                     <<"system="<<argv[arg]<<" sampled_shells="<<chosen.size()
                     <<" canonical_shell_pairs="<<pairs.size()<<" evaluated_shell_quartets="<<quartets
                     <<" ao_elements="<<ao_elements<<" spatial_tighter="<<improved
                     <<" tighter_fraction="<<(ao_elements?double(improved)/ao_elements:0.0)
                     <<" max_scaled_violation="<<max_scaled_violation
                     <<" max_cumulative_hybrid_to_cs="<<max_cumulative_ratio
                     <<" failures="<<failures<<" status="<<(failures?"FAIL":"OK")<<'\n';
            for (std::size_t t=0;t<3;++t)
                std::cout<<" threshold="<<thresholds[t]<<" cs_skips="<<cs_skips[t]
                         <<" hybrid_skips="<<hybrid_skips[t]<<'\n';
            if (failures) return 1;
        }
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
