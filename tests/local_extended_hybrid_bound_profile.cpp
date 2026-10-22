#include "chem/kumo/local_quartet_sampling.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/fchk.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    using namespace chem::kumo;
    using Clock=std::chrono::steady_clock;
    try {
        std::size_t target=1000;
        int arg=1;
        if(argc>3 && std::string(argv[1])=="--target") {
            target=std::stoull(argv[2]); arg=3;
        }
        if(arg==argc || target==0 || target>10000 || target%4)
            throw std::invalid_argument("usage: profile [--target 1000|5000|10000] file.fchk[.gz] ...");
        std::cout<<std::scientific<<std::setprecision(10);
        for(;arg<argc;++arg) {
            const std::string path=argv[arg];
            const BasisSet basis{FormattedCheckpoint(path)};
            const auto shells=local_shell_components(basis);
            const auto t0=Clock::now();
            const auto sampled=sampling::sample_quartets(shells);
            const auto t1=Clock::now();
            const auto norms=local_ao_product_norms(basis);
            const auto t2=Clock::now();
            std::vector<Eigen::Index> offsets(shells.size()+1,0);
            for(std::size_t i=0;i<shells.size();++i)
                offsets[i+1]=offsets[i]+static_cast<Eigen::Index>(shells[i].size());
            std::map<std::size_t,LocalShellPairEnvelope> cache;
            auto envelope=[&](std::size_t p)->const LocalShellPairEnvelope& {
                auto it=cache.find(p);
                if(it==cache.end()) {
                    const auto& sp=sampled.pairs.at(p);
                    it=cache.emplace(p,prepare_local_shell_pair_envelope(shells[sp.a],shells[sp.b],0.7)).first;
                }
                return it->second;
            };
            const auto t3=Clock::now();
            for(unsigned bin=0;bin<4;++bin) {
                const auto count=std::min(target/4,sampled.samples[bin].size());
                std::size_t tighter=0,invalid=0,cs_skips=0,hybrid_skips=0;
                long double work=0, saved_cs=0,saved_hybrid=0,ratio_log=0;
                std::size_t ratios=0;
                const auto b0=Clock::now();
                for(std::size_t k=0;k<count;++k) {
                    const auto& s=sampled.samples[bin][k];
                    const auto [i,j]=sampling::decode(s.id);
                    const auto& a=sampled.pairs.at(i);
                    const auto& b=sampled.pairs.at(j);
                    // Representative AO component (0,0) of each shell pair.
                    // This is NOT a bound for the entire shell quartet.
                    const auto mu=offsets[a.a], nu=offsets[a.b];
                    const auto la=offsets[b.a], si=offsets[b.b];
                    const double cs=norms(mu,nu)*norms(la,si);
                    const double spatial=local_shell_pair_envelope_bound(
                        envelope(i),0,0,envelope(j),0,0);
                    if(!std::isfinite(cs)||!std::isfinite(spatial)||cs<0||spatial<0) {
                        ++invalid;continue;
                    }
                    const double hybrid=std::min(cs,spatial);
                    tighter+=(spatial<cs);
                    if(cs>0 && hybrid>0) { ratio_log+=std::log(hybrid/cs);++ratios; }
                    work+=s.cost;
                    constexpr double threshold=1e-8;
                    if(cs<threshold) {++cs_skips;saved_cs+=s.cost;}
                    if(hybrid<threshold) {++hybrid_skips;saved_hybrid+=s.cost;}
                }
                const auto b1=Clock::now();
                std::cout<<"system="<<path<<" target="<<target<<" bin="<<bin
                         <<" sampled="<<count<<" tighter="<<tighter
                         <<" tighter_fraction="<<(count?double(tighter)/count:0.0)
                         <<" geometric_ratio_positive="<<(ratios?std::exp(ratio_log/ratios):1.0L)
                         <<" positive_ratios="<<ratios<<" invalid="<<invalid
                         <<" cs_skips="<<cs_skips<<" hybrid_skips="<<hybrid_skips
                         <<" work_proxy="<<work<<" cs_saved_work_proxy="<<saved_cs
                         <<" hybrid_saved_work_proxy="<<saved_hybrid
                         <<" bound_seconds="<<std::chrono::duration<double>(b1-b0).count()<<'\n';
                if(invalid) throw std::runtime_error("invalid bound");
            }
            const auto t4=Clock::now();
            std::cout<<"system="<<path<<" sample_seconds="<<std::chrono::duration<double>(t1-t0).count()
                     <<" norms_seconds="<<std::chrono::duration<double>(t2-t1).count()
                     <<" envelope_cache_entries="<<cache.size()
                     <<" profile_seconds="<<std::chrono::duration<double>(t4-t3).count()
                     <<" mode=representative_AO_only status=OK\n";
        }
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
