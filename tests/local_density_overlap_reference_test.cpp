// Stage 03C.8d: independent numerical validation of density-weighted AO overlap.
// Diagnostic only: the quadrature has finite accuracy and is not a production grid.
#include "chem/kumo/local_density_overlap.hpp"
#include "chem/kumo/weighted_overlap.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "lebedev_reference_grid.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
struct Rule { Eigen::VectorXd x, w; };
Rule legendre(int n) {
    Eigen::MatrixXd j=Eigen::MatrixXd::Zero(n,n);
    for(int k=1;k<n;++k) {
        const double b=double(k)/std::sqrt(4.0*k*k-1.0);
        j(k-1,k)=j(k,k-1)=b;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(j);
    if(eig.info()!=Eigen::Success) throw std::runtime_error("Gauss-Legendre diagonalization failed");
    return {(eig.eigenvalues().array()+1.0).matrix()/2.0,
            eig.eigenvectors().row(0).array().square().transpose().matrix()};
}
struct Center { std::size_t id; Eigen::Vector3d r; };
std::vector<Center> centers_of(const chem::kumo::BasisSet& basis) {
    std::vector<Center> centers;
    for(const auto& shell:basis.shells()) {
        auto it=std::find_if(centers.begin(),centers.end(),
            [&](const Center& c){return c.id==static_cast<std::size_t>(shell.center);});
        if(it==centers.end()) centers.push_back({static_cast<std::size_t>(shell.center),shell.origin});
        else if((it->r-shell.origin).norm()>1e-10)
            throw std::runtime_error("inconsistent shell centers");
    }
    return centers;
}
// Original equal-radius Becke partition: 3 applications of p(mu)=(3mu-mu^3)/2.
// No atomic-radius correction; appropriate as a controlled first diagnostic.
Eigen::VectorXd becke_weights(const Eigen::Vector3d& point,
                               const std::vector<Center>& centers) {
    const Eigen::Index n=static_cast<Eigen::Index>(centers.size());
    Eigen::VectorXd logp=Eigen::VectorXd::Zero(n);
    for(Eigen::Index a=0;a<n;++a) for(Eigen::Index b=a+1;b<n;++b) {
        const double rab=(centers[a].r-centers[b].r).norm();
        if(!(rab>1e-10)) throw std::runtime_error("coincident Becke centers");
        const double da=(point-centers[a].r).norm();
        const double db=(point-centers[b].r).norm();
        double mu=std::clamp((da-db)/rab,-1.0,1.0);
        for(int k=0;k<3;++k) mu=0.5*mu*(3.0-mu*mu);
        const double sa=std::clamp(0.5*(1.0-mu),0.0,1.0);
        const double sb=std::clamp(0.5*(1.0+mu),0.0,1.0);
        logp[a]+=sa>0.0?std::log(sa):-1e300;
        logp[b]+=sb>0.0?std::log(sb):-1e300;
    }
    const double maximum=logp.maxCoeff();
    Eigen::VectorXd p=(logp.array()-maximum).exp();
    p/=p.sum();
    return p;
}
// Angular Lebedev rules are pretabulated and independently checked below.

struct NumericalResult {
    Eigen::MatrixXd overlap;
    double rho2 = 0.0;
    std::size_t points = 0;
};

NumericalResult integrate_density(const chem::kumo::BasisEvaluator& ao,
                                  const std::vector<Center>& centers,
                                  const Eigen::MatrixXd& density,
                                  int radial_order, int angular_points) {
    const auto radial = legendre(radial_order);
    const auto angular = lebedev_reference::grid(angular_points);
    chem::kumo::WeightedOverlapAccumulator accumulator(ao);
    const Eigen::Index npoints = static_cast<Eigen::Index>(angular.size());
    Eigen::Matrix<double,3,Eigen::Dynamic> points(3,npoints);
    Eigen::VectorXd weights(npoints), rho(npoints);
    Eigen::MatrixXd values(static_cast<Eigen::Index>(ao.basis_size()),npoints);
    std::vector<double> edges{0.0,1e-5};
    for (int i=1;i<=16;++i)
        edges.push_back(1e-5*std::pow(4.0e6,double(i)/16.0));
    NumericalResult result;
    for (std::size_t a=0;a<centers.size();++a)
        for (std::size_t panel=0;panel+1<edges.size();++panel) {
            const double lo=edges[panel], hi=edges[panel+1];
            for (int ir=0;ir<radial_order;++ir) {
                const double r=lo+(hi-lo)*radial.x[ir];
                const double wr=(hi-lo)*radial.w[ir]*r*r;
                for (Eigen::Index i=0;i<npoints;++i) {
                    const auto& v=angular[static_cast<std::size_t>(i)];
                    const Eigen::Vector3d p=centers[a].r+r*Eigen::Vector3d(v.x,v.y,v.z);
                    points.col(i)=p;
                    weights[i]=wr*v.w*becke_weights(p,centers)[static_cast<Eigen::Index>(a)];
                }
                ao.evaluate_batch(points,values);
                const Eigen::MatrixXd pv = density*values;
                rho=(values.array()*pv.array()).colwise().sum().transpose().matrix();
                if (!rho.allFinite()) throw std::runtime_error("nonfinite density samples");
                result.rho2 += (weights.array()*rho.array().square()).sum();
                accumulator.add(points,weights,rho);
                result.points += static_cast<std::size_t>(npoints);
            }
        }
    result.overlap=accumulator.matrix();
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    try {
        const auto molecule=chem::kumo::Molecule::from_fchk(argv[1]);
        chem::kumo::FormattedCheckpoint fchk(argv[1]);
        chem::kumo::BasisEvaluator ao{chem::kumo::BasisSet(fchk)};
        const auto centers=centers_of(ao.basis());
        const Eigen::MatrixXd p=molecule.density_matrix(chem::kumo::Spin::Total);
        const Eigen::MatrixXd analytic=chem::kumo::local_density_overlap(ao.basis(),p);
        const double analytic_rho2=(p.cwiseProduct(analytic)).sum();
        std::cout<<std::scientific<<std::setprecision(12)
                 <<"system="<<argv[1]<<" nAO="<<ao.basis_size()
                 <<" analytic_rho2="<<analytic_rho2<<'\n';
        // Use the same radial rule and Becke partition as the previous
        // multicenter AO-overlap diagnostics; only the field changes to rho.
        for (int angular : {434,770}) {
            const auto numerical=integrate_density(ao,centers,p,20,angular);
            const double max_abs=(analytic-numerical.overlap).cwiseAbs().maxCoeff();
            const double rho2_error=std::abs(analytic_rho2-numerical.rho2);
            const double sym=(numerical.overlap-numerical.overlap.transpose()).cwiseAbs().maxCoeff();
            std::cout<<"lebedev="<<angular<<" points="<<numerical.points
                     <<" max_abs="<<max_abs
                     <<" rho2_quadrature="<<numerical.rho2
                     <<" rho2_abs_error="<<rho2_error
                     <<" symmetry="<<sym<<'\n';
            if (!std::isfinite(max_abs) || !std::isfinite(rho2_error) ||
                !std::isfinite(sym) || sym>1e-9)
                throw std::runtime_error("invalid numerical density overlap");
        }
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
    return 0;
}
