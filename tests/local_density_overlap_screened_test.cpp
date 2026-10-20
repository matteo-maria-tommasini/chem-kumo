#include "chem/kumo/local_density_overlap_screened.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    using namespace chem::kumo;
    if (argc<2) return 2;
    try {
        for (int arg=1;arg<argc;++arg) {
            const auto molecule=Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const Eigen::MatrixXd density=molecule.density_matrix(Spin::Total);
            const auto t0=std::chrono::steady_clock::now();
            const Eigen::MatrixXd reference=local_density_overlap_symmetric_paired(basis,density);
            const double reference_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
            for (const double eps:{0.0,1e-6,1e-8,1e-10}) {
                const auto result=local_density_overlap_screened(basis,density,eps);
                const auto& stats=result.statistics;
                const double err=(reference-result.matrix).cwiseAbs().maxCoeff();
                const double sym=(result.matrix-result.matrix.transpose()).cwiseAbs().maxCoeff();
                const double total=stats.norm_seconds+stats.selection_seconds+stats.pair_preparation_seconds+stats.contraction_seconds;
                std::cout<<std::scientific<<std::setprecision(12)
                         <<"system="<<argv[arg]<<" epsilon="<<eps
                         <<" quartets="<<stats.total_quartets
                         <<" skipped="<<stats.skipped_quartets
                         <<" max_bound="<<stats.max_error_bound
                         <<" max_abs="<<err<<" symmetry="<<sym
                         <<" reference_seconds="<<reference_seconds
                         <<" norm_seconds="<<stats.norm_seconds
                         <<" selection_seconds="<<stats.selection_seconds
                         <<" preparation_seconds="<<stats.pair_preparation_seconds
                         <<" contraction_seconds="<<stats.contraction_seconds
                         <<" total_seconds="<<total
                         <<" speedup="<<reference_seconds/total<<'\n';
                if (!result.matrix.allFinite() || sym>1e-12 ||
                    err>eps+2e-11*std::max(1.0,reference.cwiseAbs().maxCoeff()) ||
                    stats.max_error_bound>eps*(1.0+1e-12))
                    throw std::runtime_error("screened density overlap failed");
            }
        }
    } catch (const std::exception& e) {
        std::cerr<<"screened density overlap: "<<e.what()<<'\n';
        return 1;
    }
}
