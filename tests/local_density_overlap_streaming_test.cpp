#include "chem/kumo/local_density_overlap_streaming.hpp"
#include "chem/kumo/local_density_overlap_symmetric_paired.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    if (argc<2) return 2;
    try {
        using namespace chem::kumo;
        std::cout<<std::scientific<<std::setprecision(12);
        for (int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule=Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd density=molecule.density_matrix(Spin::Total);
            const Eigen::MatrixXd reference=local_density_overlap_symmetric_paired(basis,density);
            for (double eps : {0.0,1e-8,1e-10}) {
                const auto result=local_density_overlap_streaming(basis,density,eps,1024*1024);
                const double err=(result.matrix-reference).cwiseAbs().maxCoeff();
                const double tol=eps+2e-11*std::max(1.0,reference.cwiseAbs().maxCoeff());
                if (err>tol || result.statistics.max_error_bound>eps+1e-18)
                    throw std::runtime_error("streaming density overlap mismatch");
                std::cout<<"system="<<argv[arg]<<" epsilon="<<eps
                         <<" max_abs="<<err
                         <<" max_bound="<<result.statistics.max_error_bound
                         <<" skipped="<<result.statistics.skipped_quartets
                         <<" pair_preparations="<<result.statistics.prepared_pair_instances
                         <<" peak_pair_bytes="<<result.statistics.peak_pair_bytes
                         <<" contraction_seconds="<<result.statistics.contraction_seconds<<'\n';
            }
        }
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
