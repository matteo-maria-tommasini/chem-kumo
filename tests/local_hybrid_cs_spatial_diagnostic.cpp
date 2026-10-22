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
    if (argc < 2) {
        std::cerr << "usage: diagnostic file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const std::size_t n = std::min<std::size_t>(6, shells.size());
            std::vector<Eigen::Index> offsets(shells.size()+1, 0);
            for (std::size_t a=0; a<shells.size(); ++a)
                offsets[a+1]=offsets[a]+static_cast<Eigen::Index>(shells[a].size());
            const auto norms = local_ao_product_norms(basis);
            std::vector<LocalShellPairEnvelope> env;
            env.reserve(n*n);
            for (std::size_t a=0; a<n; ++a)
                for (std::size_t b=0; b<n; ++b)
                    env.push_back(prepare_local_shell_pair_envelope(shells[a],shells[b],0.7));
            constexpr std::array<double,4> thresholds{1e-6,1e-8,1e-10,1e-12};
            std::array<std::size_t,4> cs_skips{}, hybrid_skips{};
            std::size_t quartets=0, elements=0, spatial_tighter=0, weighted_checks=0;
            long double max_excess=0, max_scaled_excess=0;
            long double log_ratio=0;
            std::size_t ratio_count=0;
            auto start=std::chrono::steady_clock::now();
            // Ordered shell quartets deliberately exercise both contraction directions.
            for (std::size_t a=0; a<n; ++a)
                for (std::size_t b=0; b<n; ++b)
                    for (std::size_t c=0; c<n; ++c)
                        for (std::size_t d=0; d<n; ++d) {
                            const auto& ab=env[a*n+b];
                            const auto& cd=env[c*n+d];
                            const auto q=local_shell_quartet({&shells[a],&shells[b],&shells[c],&shells[d]});
                            ++quartets;
                            for (std::size_t i=0;i<ab.na;++i)
                                for (std::size_t j=0;j<ab.nb;++j)
                                    for (std::size_t k=0;k<cd.na;++k)
                                        for (std::size_t l=0;l<cd.nb;++l) {
                                            const auto mu=offsets[a]+static_cast<Eigen::Index>(i);
                                            const auto nu=offsets[b]+static_cast<Eigen::Index>(j);
                                            const auto lambda=offsets[c]+static_cast<Eigen::Index>(k);
                                            const auto sigma=offsets[d]+static_cast<Eigen::Index>(l);
                                            const double cs=norms(mu,nu)*norms(lambda,sigma);
                                            const double spatial=local_shell_pair_envelope_bound(ab,i,j,cd,k,l);
                                            if (!(cs>=0) || !(spatial>=0) || !std::isfinite(cs) || !std::isfinite(spatial))
                                                throw std::runtime_error("nonfinite or negative bound");
                                            const double hybrid=std::min(cs,spatial);
                                            const double exact=std::abs(q.at(i,j,k,l));
                                            const long double excess=std::max(0.0L,static_cast<long double>(exact)-hybrid);
                                            const long double scale=std::max({1.0L,static_cast<long double>(exact),static_cast<long double>(hybrid)});
                                            max_excess=std::max(max_excess,excess);
                                            max_scaled_excess=std::max(max_scaled_excess,excess/scale);
                                            if (excess>2e-11L*scale)
                                                throw std::runtime_error("hybrid bound violation");
                                            spatial_tighter+=(spatial<cs);
                                            if (cs>0 && hybrid>0) {
                                                log_ratio+=std::log(static_cast<long double>(hybrid)/cs);
                                                ++ratio_count;
                                            }
                                            for (std::size_t t=0;t<thresholds.size();++t) {
                                                cs_skips[t]+=(cs<thresholds[t]);
                                                hybrid_skips[t]+=(hybrid<thresholds[t]);
                                            }
                                            // Deterministic symmetric nonnegative AO weights. Both
                                            // canonical pair multiplicities are checked independently.
                                            const long double w1=1.0L/(1.0L+std::abs(lambda-sigma));
                                            const long double w2=1.0L/(1.0L+std::abs(mu-nu));
                                            const long double m1=(lambda==sigma ? 1.0L : 2.0L);
                                            const long double m2=(mu==nu ? 1.0L : 2.0L);
                                            for (const auto weight : {w1*m1,w2*m2}) {
                                                if (weight*exact>weight*hybrid+2e-11L*std::max(1.0L,weight*exact))
                                                    throw std::runtime_error("weighted bound violation");
                                                ++weighted_checks;
                                            }
                                            ++elements;
                                        }
                        }
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            std::cout<<std::scientific<<std::setprecision(12)
                     <<"system="<<argv[arg]<<" sampled_shells="<<n
                     <<" ordered_quartets="<<quartets<<" ao_elements="<<elements
                     <<" weighted_checks="<<weighted_checks
                     <<" spatial_tighter="<<spatial_tighter
                     <<" tighter_fraction="<<(elements?double(spatial_tighter)/elements:0)
                     <<" geometric_mean_hybrid_to_cs="<<(ratio_count?std::exp(log_ratio/ratio_count):1.0L)
                     <<" max_positive_violation="<<max_excess
                     <<" max_scaled_violation="<<max_scaled_excess
                     <<" seconds="<<seconds<<" status=OK\n";
            for (std::size_t t=0;t<thresholds.size();++t)
                std::cout<<" threshold="<<thresholds[t]<<" cs_skips="<<cs_skips[t]
                         <<" hybrid_skips="<<hybrid_skips[t]<<'\n';
        }
    } catch (const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
