#include "chem/kumo/local_density_overlap_plan.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

namespace {
std::uint64_t resident_bytes() {
#ifdef __APPLE__
    mach_task_basic_info_data_t info{};
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
        return static_cast<std::uint64_t>(info.resident_size);
#endif
    return 0;
}
}

int main(int argc, char** argv) {
    using namespace chem::kumo;
    if (argc < 2) {
        std::cerr << "usage: memory-test [--build-max-mib N] file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        double max_mib = 0.0; // Estimate-only by default, safe for HBC.
        std::vector<std::string> paths;
        for (int i=1; i<argc; ++i) {
            const std::string arg=argv[i];
            if (arg=="--build-max-mib") {
                if (++i==argc) throw std::invalid_argument("missing MiB limit");
                max_mib=std::stod(argv[i]);
                if (!(max_mib>0.0) || !std::isfinite(max_mib))
                    throw std::invalid_argument("invalid MiB limit");
            } else paths.push_back(arg);
        }
        for (const auto& path:paths) {
            const auto t0=std::chrono::steady_clock::now();
            const BasisSet basis{FormattedCheckpoint(path)};
            const auto shells=local_shell_components(basis);
            std::size_t nao=0, total_primitives=0;
            std::vector<std::size_t> component_primitive_counts;
            for (const auto& shell:shells) {
                nao+=shell.size();
                for (const auto& component:shell) {
                    component_primitive_counts.push_back(component.size());
                    total_primitives+=component.size();
                }
            }
            // Canonical shell pairs contain all ordered AO-component products
            // within diagonal shell blocks, as in prepare_local_shell_pair().
            std::uint64_t pair_count=0, component_vectors=0, max_pair=0;
            std::size_t offset_a=0;
            for (std::size_t a=0;a<shells.size();++a) {
                std::size_t offset_b=offset_a;
                for (std::size_t b=a;b<shells.size();++b) {
                    std::uint64_t local=0;
                    for (std::size_t i=0;i<shells[a].size();++i)
                        for (std::size_t j=0;j<shells[b].size();++j) {
                            const auto ni=component_primitive_counts[offset_a+i];
                            const auto nj=component_primitive_counts[offset_b+j];
                            local+=static_cast<std::uint64_t>(ni)*nj;
                            ++component_vectors;
                        }
                    pair_count+=local;
                    max_pair=std::max(max_pair,local);
                    offset_b+=shells[b].size();
                }
                offset_a+=shells[a].size();
            }
            const auto shell_pairs=static_cast<std::uint64_t>(shells.size())*(shells.size()+1)/2;
            const auto quartets=shell_pairs*(shell_pairs+1)/2;
            const auto payload=pair_count*sizeof(GaussianPairIntermediate);
            const auto norms=static_cast<std::uint64_t>(nao)*nao*sizeof(double);
            const auto vectors=component_vectors*sizeof(std::vector<GaussianPairIntermediate>);
            const auto pair_objects=shell_pairs*sizeof(LocalPreparedShellPair);
            const auto pair_meta=shell_pairs*(2*sizeof(std::size_t)+sizeof(double));
            const auto lower=payload+norms+vectors+pair_objects+pair_meta;
            const double mib=1048576.0;
            std::cout<<std::setprecision(10)<<"system="<<path
                     <<" shells="<<shells.size()<<" nao="<<nao
                     <<" shell_pairs="<<shell_pairs<<" quartets="<<quartets
                     <<" component_primitive_entries="<<total_primitives
                     <<" primitive_pairs="<<pair_count<<" max_shell_pair_primitives="<<max_pair
                     <<" primitive_payload_mib="<<payload/mib
                     <<" norms_mib="<<norms/mib
                     <<" component_vector_headers_mib="<<vectors/mib
                     <<" estimated_lower_mib="<<lower/mib
                     <<" estimate_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
            if (max_mib>0.0 && lower/mib<=max_mib) {
                const auto before=resident_bytes();
                const auto t1=std::chrono::steady_clock::now();
                {
                    const LocalDensityOverlapPlan plan(basis);
                    const auto after=resident_bytes();
                    std::cout<<" plan_primitive_pairs="<<plan.primitive_pair_count()
                             <<" setup_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count()
                             <<" resident_before_mib="<<before/mib
                             <<" resident_with_plan_mib="<<after/mib;
                    if (plan.primitive_pair_count()!=pair_count)
                        throw std::runtime_error("estimated primitive pair count differs from plan");
                }
            } else std::cout<<" plan_build=skipped";
            std::cout<<'\n';
        }
    } catch (const std::exception& e) {
        std::cerr<<"memory diagnostic: "<<e.what()<<'\n';
        return 1;
    }
}
