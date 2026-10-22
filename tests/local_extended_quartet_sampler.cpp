#include "chem/kumo/local_quartet_sampling.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        if(argc<2) throw std::invalid_argument("usage: kumo-local-extended-quartet-sampler file.fchk[.gz] ...");
        std::cout<<std::scientific<<std::setprecision(10);
        for(int arg=1;arg<argc;++arg) {
            const std::string path=argv[arg];
            const chem::kumo::BasisSet basis{chem::kumo::FormattedCheckpoint(path)};
            const auto shells=chem::kumo::local_shell_components(basis);
            const std::uint64_t ns=shells.size(), np=chem::kumo::sampling::triangular(ns), nq=chem::kumo::sampling::triangular(np);
            const auto result=chem::kumo::sampling::sample_quartets(shells);
            const auto& samples=result.samples;
            const auto unique_draws=result.unique_draws;
            std::cout<<"system="<<path<<" shells="<<ns<<" shell_pairs="<<np
                     <<" formal_quartets="<<nq<<" unique_draws="<<unique_draws<<'\n';
            for(const std::size_t target:{1000,5000,10000}) {
                const std::size_t quota=target/4;
                std::size_t selected=0;
                long double total_work=0;
                for(unsigned bin=0;bin<4;++bin) {
                    const auto count=std::min(quota,samples[bin].size());
                    long double work=0;
                    double min_d=std::numeric_limits<double>::infinity(),max_d=0;
                    for(std::size_t k=0;k<count;++k) {
                        const auto& s=samples[bin][k];
                        work+=s.cost; min_d=std::min(min_d,s.distance); max_d=std::max(max_d,s.distance);
                    }
                    selected+=count;total_work+=work;
                    std::cout<<" target="<<target<<" bin="<<bin<<" sampled="<<count
                             <<" requested="<<quota<<" work_proxy="<<static_cast<double>(work)
                             <<" min_distance="<<(count?min_d:0)
                             <<" max_distance="<<max_d<<'\n';
                }
                std::cout<<" target="<<target<<" selected="<<selected
                         <<" underfilled="<<(target-selected)
                         <<" total_work_proxy="<<static_cast<double>(total_work)<<'\n';
            }
            std::cout<<"system="<<path<<" status=OK\n";
        }
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
