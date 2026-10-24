#include "chem/kumo/local_density_pair_graph.hpp"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        using namespace chem::kumo;
        // Real molecular AO order and overlap, no four-center integrals.
        for(int k=1;k<argc;++k) {
            const FormattedCheckpoint fchk(argv[k]);
            const BasisSet basis(fchk);
            const Molecule mol=Molecule::from_fchk(argv[k]);
            const Eigen::MatrixXd p=mol.density_matrix(Spin::Total);
            const Eigen::MatrixXd s=BasisEvaluator(basis).overlap_matrix();
            const auto scores=local_density_pair_graph(basis,mol,p,s);
            const auto n=mol.atoms().size();
            if(scores.size()!=n*(n-1)/2) throw std::runtime_error("pair count mismatch");
            for(const auto& x:scores) {
                if(!std::isfinite(x.distance_bohr)||!std::isfinite(x.raw_frobenius)||
                   !std::isfinite(x.metric_frobenius)||!std::isfinite(x.normalized)||
                   x.raw_frobenius<0||x.metric_frobenius<0||x.normalized<0)
                    throw std::runtime_error("invalid score");
                std::cout<<"pair,"<<argv[k]<<','<<x.a<<','<<x.b<<','<<x.distance_bohr<<','
                         <<x.raw_frobenius<<','<<x.metric_frobenius<<','<<x.normalized<<'\n';
            }
            std::cout<<"system="<<argv[k]<<" atoms="<<n<<" pairs="<<scores.size()<<" status=OK\n";
        }
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
