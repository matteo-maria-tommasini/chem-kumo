#include "chem/kumo/local_density_overlap_plan.hpp"
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
            const auto total=molecule.density_matrix(Spin::Total);
            const auto spin=molecule.density_matrix(Spin::Difference);
            const auto t0=std::chrono::steady_clock::now();
            const LocalDensityOverlapPlan plan(basis);
            const double setup=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-t0).count();
            for (const double eps:{0.0,1e-8,1e-10}) {
                for (int channel=0;channel<2;++channel) {
                    const Eigen::MatrixXd& p=channel==0?total:spin;
                    const auto baseline=local_density_overlap_screened(basis,p,eps);
                    const auto t1=std::chrono::steady_clock::now();
                    const auto reused=plan.evaluate(p,eps);
                    const double seconds=std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-t1).count();
                    const double err=(baseline.matrix-reused.matrix).cwiseAbs().maxCoeff();
                    const double scale=std::max(1.0,baseline.matrix.cwiseAbs().maxCoeff());
                    std::cout<<std::scientific<<std::setprecision(12)
                             <<"system="<<argv[arg]<<" channel="<<(channel==0?"total":"spin")
                             <<" epsilon="<<eps<<" setup_seconds="<<setup
                             <<" evaluation_seconds="<<seconds
                             <<" primitive_pairs="<<plan.primitive_pair_count()
                             <<" skipped="<<reused.statistics.skipped_quartets
                             <<" max_bound="<<reused.statistics.max_error_bound
                             <<" max_abs="<<err<<'\n';
                    if (!reused.matrix.allFinite() || err>2e-11*scale ||
                        reused.statistics.skipped_quartets!=baseline.statistics.skipped_quartets ||
                        reused.statistics.max_error_bound>eps*(1.0+1e-12))
                        throw std::runtime_error("prepared plan mismatch");
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr<<"prepared plan: "<<e.what()<<'\n';
        return 1;
    }
}
