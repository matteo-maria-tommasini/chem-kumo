#include "chem/kumo/weighted_overlap.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
namespace {
struct Rule {
    Eigen::VectorXd nodes, weights;
    explicit Rule(int n) {
        Eigen::MatrixXd j=Eigen::MatrixXd::Zero(n,n);
        for(int k=1;k<n;++k) j(k-1,k)=j(k,k-1)=std::sqrt(k/2.0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(j);
        if(eig.info()!=Eigen::Success) throw std::runtime_error("Gauss-Hermite eigensolver failed");
        nodes=eig.eigenvalues();
        weights=(eig.eigenvectors().row(0).array().square()*std::sqrt(std::acos(-1.0))).matrix().transpose();
    }
};
// Diagnostic-only global tensor Gauss-Hermite quadrature. This is NOT a
// general molecular quadrature: the Gaussian envelope is user-selected.
// The 3D nodes are generated and accumulated in slices to bound memory.
double integrate_ao_pair(const BasisEvaluator& ao, Eigen::Index mu, Eigen::Index nu,
                         const Eigen::Vector3d& center, double exponent, int order) {
    Rule rule(order);
    WeightedOverlapAccumulator accum(ao);
    Eigen::Matrix<double,3,Eigen::Dynamic> xyz(3,order);
    Eigen::VectorXd q(order), field=Eigen::VectorXd::Ones(order);
    const double inv=1.0/std::sqrt(exponent);
    for(int i=0;i<order;++i) for(int j=0;j<order;++j) {
        for(int k=0;k<order;++k) {
            const Eigen::Vector3d t(rule.nodes[i],rule.nodes[j],rule.nodes[k]);
            xyz.col(k)=center+inv*t;
            // Undo the Hermite envelope to integrate over Lebesgue measure.
            q[k]=rule.weights[i]*rule.weights[j]*rule.weights[k]
                 *std::exp(t.squaredNorm())*inv*inv*inv;
        }
        accum.add(xyz,q,field);
    }
    return accum.matrix()(mu,nu);
}
}
int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        for(int a=1;a<argc;++a) {
            FormattedCheckpoint fchk(argv[a]);
            BasisEvaluator ao{BasisSet(fchk)};
            const auto& shell=ao.basis().shells().front();
            if(shell.exponents.empty()) throw std::runtime_error("empty first shell");
            const auto minmax=std::minmax_element(shell.exponents.begin(),shell.exponents.end());
            // Geometric mean balances tight and diffuse primitives in the first shell.
            const double exponent=2.0*std::sqrt((*minmax.first)*(*minmax.second));
            if(!(exponent>0.0) || !std::isfinite(exponent)) throw std::runtime_error("invalid envelope");
            const Eigen::MatrixXd exact=ao.overlap_matrix();
            const Eigen::Index mu=0,nu=0;
            std::cout<<std::setprecision(12)<<argv[a]<<" nAO="<<ao.basis_size()
                     <<" envelope="<<exponent<<" reference="<<exact(mu,nu);
            for(int order : {8,12,16,24}) {
                const double value=integrate_ao_pair(ao,mu,nu,shell.origin,exponent,order);
                const double error=std::abs(value-exact(mu,nu));
                if(!std::isfinite(value) || !std::isfinite(error))
                    throw std::runtime_error("non-finite quadrature result");
                std::cout<<" n="<<order<<" error="<<error;
            }
            std::cout<<'\n';
        }
    } catch(const std::exception& e) {
        std::cerr<<"weighted-overlap quadrature diagnostic: "<<e.what()<<'\n';
        return 1;
    }
}
