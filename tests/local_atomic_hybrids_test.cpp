#include "chem/kumo/local_atomic_hybrids.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc,char** argv) {
    try {
        using namespace chem::kumo;
        // Synthetic rank-deficient metric: check rank truncation and metric normalization.
        LocalAtomicW synthetic;
        synthetic.ao_indices={0,1,2};
        Eigen::MatrixXd s=Eigen::MatrixXd::Identity(3,3);
        s(2,2)=1e-14;
        Eigen::MatrixXd w=Eigen::MatrixXd::Zero(3,3);
        w(0,0)=2.; w(1,1)=1.; w(2,2)=1e-14;
        synthetic.matrices={w};
        const auto reduced=local_atomic_hybrids(synthetic,s,1e-10);
        if(reduced.rank!=2 || reduced.channels[0].orthogonality_max>1e-10)
            throw std::runtime_error("rank truncation regression");
        if(argc<2) return 2;
        for(int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule=Molecule::from_fchk(argv[arg]);
            const std::vector<Eigen::MatrixXd> densities{
                molecule.density_matrix(Spin::Total),molecule.density_matrix(Spin::Difference)};
            const Eigen::MatrixXd overlap=BasisEvaluator(basis).overlap_matrix();
            const std::size_t center=basis.shells().front().center;
            const auto block=local_atomic_w(basis,densities,center);
            const auto result=local_atomic_hybrids(block,overlap);
            if(result.channels.size()!=2 || result.rank<1)
                throw std::runtime_error("missing hybrid channels");
            for(std::size_t channel=0;channel<result.channels.size();++channel) {
                const auto& ch=result.channels[channel];
                if(ch.orthogonality_max>1e-8 || ch.residual_max>1e-8)
                    throw std::runtime_error("generalized eigensystem regression");
                for(Eigen::Index i=0;i<result.rank;++i) {
                    const double expectation=(ch.coefficients.col(i).transpose()*block.matrices[channel]*ch.coefficients.col(i))(0,0);
                    if(std::abs(expectation-ch.eigenvalues[i])>1e-8*std::max(1.,std::abs(expectation)))
                        throw std::runtime_error("Rayleigh quotient mismatch");
                }
            }
            std::cout<<"system="<<argv[arg]<<" center="<<center<<" nao="<<block.ao_indices.size()
                     <<" rank="<<result.rank<<" eigenvalues_total="<<result.channels[0].eigenvalues.transpose()
                     <<" status=OK\n";
        }
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
