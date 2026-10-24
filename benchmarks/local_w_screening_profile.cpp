#include "chem/kumo/local_density_overlap_multichannel_hierarchical.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace chem::kumo;
int main(int argc, char** argv) {
    try {
        double epsilon = 1e-8;
        std::vector<std::string> files;
        for (int i=1; i<argc; ++i) {
            std::string a=argv[i];
            if (a=="--epsilon" && i+1<argc) epsilon=std::stod(argv[++i]);
            else if (a=="--help") {
                std::cerr << "usage: kumo-local-w-screening-profile [--epsilon E] file.fchk[.gz] ...\n";
                return 0;
            } else if (a.rfind("--",0)==0) throw std::invalid_argument("unknown option: "+a);
            else files.push_back(a);
        }
        if (files.empty() || !std::isfinite(epsilon) || epsilon<0)
            throw std::invalid_argument("invalid arguments");
        std::cout << "system,shell_pairs,formal_quartets,selected_quartets,skipped_quartets,norm_seconds,ordering_seconds,total_seconds,max_error_bound\n";
        std::cout << std::setprecision(16);
        for (const auto& file : files) {
            const BasisSet basis{FormattedCheckpoint(file)};
            const auto molecule=Molecule::from_fchk(file);
            std::vector<Eigen::MatrixXd> densities{molecule.density_matrix(Spin::Total)};
            auto result=local_density_overlap_multichannel_hierarchical(
                basis,densities,{epsilon},64*1024*1024,
                MultichannelScreeningOrder::MaximumImportance,true);
            const auto& s=result.statistics;
            if (s.evaluated_quartets!=0 ||
                s.candidate_visits+s.skipped_quartets!=s.total_quartets)
                throw std::runtime_error("inconsistent screening-only counters");
            std::cout << file << ',' << s.shell_pairs << ',' << s.total_quartets
                << ',' << s.candidate_visits << ',' << s.skipped_quartets
                << ',' << s.norm_seconds << ',' << s.ordering_seconds
                << ',' << s.total_seconds << ',' << s.max_error_bounds.at(0) << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "screening profile error: " << e.what() << '\n'; return 1;
    }
}
