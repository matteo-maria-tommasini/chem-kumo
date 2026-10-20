#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
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

int main(int argc, char** argv) {
    using namespace chem::kumo;
    using clock = std::chrono::steady_clock;
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const std::size_t n = std::min<std::size_t>(6, shells.size());
            const auto t0 = clock::now();
            const Eigen::MatrixXd norms = local_ao_product_norms(basis);
            const auto t1 = clock::now();
            std::vector<Eigen::Index> offsets(shells.size()+1, 0);
            for (std::size_t a=0; a<shells.size(); ++a)
                offsets[a+1]=offsets[a]+static_cast<Eigen::Index>(shells[a].size());
            std::vector<LocalShellPairEnvelope> envelopes;
            envelopes.reserve(n*n);
            for (std::size_t a=0; a<n; ++a)
                for (std::size_t b=0; b<n; ++b)
                    envelopes.push_back(prepare_local_shell_pair_envelope(shells[a],shells[b],0.7));
            const auto t2 = clock::now();
            std::array<std::size_t,4> cs_skip{}, combined_skip{};
            constexpr std::array<double,4> thresholds{1e-6,1e-8,1e-10,1e-12};
            std::size_t elements=0, tighter=0, quartets=0;
            double violation=0.0, max_ratio=0.0;
            long double log_ratio_sum=0.0L;
            std::size_t ratio_count=0;
            double envelope_seconds=0.0, integral_seconds=0.0;
            volatile double sink=0.0;
            for (std::size_t a=0; a<n; ++a)
                for (std::size_t b=0; b<n; ++b)
                    for (std::size_t c=0; c<n; ++c)
                        for (std::size_t d=0; d<n; ++d) {
                            const auto& ab=envelopes[a*n+b];
                            const auto& cd=envelopes[c*n+d];
                            const auto q0=clock::now();
                            const auto q=local_shell_quartet({&shells[a],&shells[b],&shells[c],&shells[d]});
                            const auto q1=clock::now();
                            integral_seconds+=std::chrono::duration<double>(q1-q0).count();
                            ++quartets;
                            for (std::size_t i=0;i<ab.na;++i)
                                for (std::size_t j=0;j<ab.nb;++j)
                                    for (std::size_t k=0;k<cd.na;++k)
                                        for (std::size_t l=0;l<cd.nb;++l) {
                                            const double cs=norms(offsets[a]+static_cast<Eigen::Index>(i),offsets[b]+static_cast<Eigen::Index>(j))*
                                                norms(offsets[c]+static_cast<Eigen::Index>(k),offsets[d]+static_cast<Eigen::Index>(l));
                                            const auto e0=clock::now();
                                            const double env=local_shell_pair_envelope_bound(ab,i,j,cd,k,l);
                                            const auto e1=clock::now();
                                            envelope_seconds+=std::chrono::duration<double>(e1-e0).count();
                                            const double exact=std::abs(q.at(i,j,k,l));
                                            if (!std::isfinite(env) || !std::isfinite(cs)) throw std::runtime_error("nonfinite bound");
                                            violation=std::max(violation,exact-std::min(cs,env));
                                            if (exact>std::min(cs,env)+2e-11*std::max({1.0,exact,cs,env}))
                                                throw std::runtime_error("bound violation");
                                            if (env<cs) ++tighter;
                                            if (cs>0 && env>0) {
                                                const double ratio=env/cs;
                                                if (std::isfinite(ratio) && ratio>0) {
                                                    log_ratio_sum+=std::log(ratio);
                                                    ++ratio_count;
                                                    max_ratio=std::max(max_ratio,ratio);
                                                }
                                            }
                                            for (std::size_t t=0;t<thresholds.size();++t) {
                                                if (cs<thresholds[t]) ++cs_skip[t];
                                                if (std::min(cs,env)<thresholds[t]) ++combined_skip[t];
                                            }
                                            sink=sink+exact*1e-15;
                                            ++elements;
                                        }
                        }
            std::cout<<std::scientific<<std::setprecision(12)
                <<"system="<<argv[arg]<<" shells="<<n<<" quartets="<<quartets
                <<" elements="<<elements<<" envelope_tighter="<<tighter
                <<" tighter_fraction="<<static_cast<double>(tighter)/elements
                <<" geometric_mean_envelope_to_cs="<<(ratio_count?std::exp(log_ratio_sum/ratio_count):0.0L)
                <<" max_envelope_to_cs="<<max_ratio
                <<" max_positive_violation="<<violation
                <<" norms_seconds="<<std::chrono::duration<double>(t1-t0).count()
                <<" preparation_seconds="<<std::chrono::duration<double>(t2-t1).count()
                <<" integral_seconds="<<integral_seconds
                <<" envelope_seconds="<<envelope_seconds
                <<" cost_ratio="<<envelope_seconds/integral_seconds
                <<" checksum="<<sink<<'\n';
            for (std::size_t t=0;t<thresholds.size();++t)
                std::cout<<" threshold="<<thresholds[t]<<" cs_skips="<<cs_skip[t]
                         <<" combined_skips="<<combined_skip[t]<<'\n';
        }
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
