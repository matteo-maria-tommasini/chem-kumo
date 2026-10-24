#include "chem/kumo/local_density_overlap_multichannel_hierarchical.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef __APPLE__
#include <sys/resource.h>
#endif
using namespace chem::kumo;
using Clock = std::chrono::steady_clock;
static double seconds(Clock::time_point a) { return std::chrono::duration<double>(Clock::now()-a).count(); }
static double peak_rss_mb() {
#ifdef __APPLE__
    struct rusage r{};
    if (getrusage(RUSAGE_SELF, &r) == 0) return double(r.ru_maxrss)/1048576.0;
#endif
    return -1.0;
}
static void help() {
    std::cerr << "usage: kumo-local-w-backend-benchmark [--repeats N] [--epsilon E] [--memory-mib N] [--reference] file.fchk[.gz] ...\n";
}
int main(int argc,char** argv) {
    try {
        int repeats=3; double epsilon=1e-8; std::size_t memory_mib=64;
        bool reference=false; std::vector<std::string> files;
        for(int i=1;i<argc;++i) {
            const std::string a=argv[i];
            if(a=="--repeats" && i+1<argc) repeats=std::stoi(argv[++i]);
            else if(a=="--epsilon" && i+1<argc) epsilon=std::stod(argv[++i]);
            else if(a=="--memory-mib" && i+1<argc) memory_mib=std::stoull(argv[++i]);
            else if(a=="--reference") reference=true;
            else if(a=="--help") { help(); return 0; }
            else if(a.rfind("--",0)==0) throw std::invalid_argument("unknown option: "+a);
            else files.push_back(a);
        }
        if(files.empty() || repeats<1 || !std::isfinite(epsilon) || epsilon<0 || memory_mib==0 || memory_mib>std::size_t(-1)/(1024*1024)) { help(); return 2; }
        const std::size_t cap=memory_mib*1024*1024;
        std::cout << "system,backend,repetition,wall_seconds,internal_seconds,norm_seconds,ordering_seconds,contraction_seconds,total_quartets,evaluated_quartets,skipped_quartets,peak_pair_mib,process_peak_rss_mib,max_error_bound,checksum\n";
        std::cout << std::setprecision(16);
        for(const auto& file:files) {
            const BasisSet basis{FormattedCheckpoint(file)};
            const auto molecule=Molecule::from_fchk(file);
            const Eigen::MatrixXd density=molecule.density_matrix(Spin::Total);
            const std::vector<Eigen::MatrixXd> channels{density};
            const std::vector<double> budgets{epsilon};
            for(int rep=0;rep<repeats;++rep) {
                if(reference) {
                    const auto t=Clock::now();
                    const auto r=local_density_overlap_multichannel(basis,channels);
                    const double wall=seconds(t);
                    std::cout<<file<<",reference,"<<rep<<','<<wall<<','<<r.statistics.total_seconds
                      <<",0,0,"<<r.statistics.quartet_seconds<<",0,"<<r.statistics.evaluated_quartets
                      <<",0,0,"<<peak_rss_mb()<<",0,"<<r.matrices[0].sum()<<'\n';
                }
                {
                    const auto t=Clock::now();
                    const auto r=local_density_overlap_multichannel_screened(basis,channels,budgets,cap,MultichannelScreeningOrder::Canonical);
                    const double wall=seconds(t); const auto& s=r.statistics;
                    std::cout<<file<<",screened,"<<rep<<','<<wall<<','<<s.total_seconds<<','<<s.norm_seconds<<','<<s.ordering_seconds<<','<<s.contraction_seconds<<','<<s.total_quartets<<','<<s.evaluated_quartets<<','<<s.skipped_quartets<<','<<double(s.peak_pair_bytes)/1048576.0<<','<<peak_rss_mb()<<','<<s.max_error_bounds.at(0)<<','<<r.matrices[0].sum()<<'\n';
                }
                {
                    const auto t=Clock::now();
                    const auto r=local_density_overlap_multichannel_hierarchical(basis,channels,budgets,cap);
                    const double wall=seconds(t); const auto& s=r.statistics;
                    std::cout<<file<<",hierarchical,"<<rep<<','<<wall<<','<<s.total_seconds<<','<<s.norm_seconds<<','<<s.ordering_seconds<<','<<s.contraction_seconds<<','<<s.total_quartets<<','<<s.evaluated_quartets<<','<<s.skipped_quartets<<','<<double(s.peak_pair_bytes)/1048576.0<<','<<peak_rss_mb()<<','<<s.max_error_bounds.at(0)<<','<<r.matrices[0].sum()<<'\n';
                }
            }
        }
    } catch(const std::exception& e) { std::cerr<<"benchmark error: "<<e.what()<<'\n'; return 1; }
}
