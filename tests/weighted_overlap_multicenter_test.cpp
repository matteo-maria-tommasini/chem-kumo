// Stage 03C.6: diagnostic molecular Becke quadrature against analytic AO overlap.
// Test-only quadrature. No production grid policy or density backend changes.
#include "chem/kumo/weighted_overlap.hpp"
#include "chem/kumo/fchk.hpp"
#include "lebedev_reference_grid.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <utility>
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
std::vector<std::size_t> ao_centers(const chem::kumo::BasisSet& basis) {
    std::vector<std::size_t> result;
    for (const auto& shell : basis.shells())
        result.insert(result.end(), shell.function_count(), shell.center);
    if (result.size() != basis.function_count())
        throw std::runtime_error("AO center map size mismatch");
    return result;
}
// Angular Lebedev rules are pretabulated and independently checked below.
Eigen::MatrixXd integrate_lebedev(const chem::kumo::BasisEvaluator& ao,
                                 const std::vector<Center>& centers,
                                 int radial_order, int angular_points) {
    const auto radial=legendre(radial_order);
    const auto angular=lebedev_reference::grid(angular_points);
    chem::kumo::WeightedOverlapAccumulator accumulator(ao);
    const Eigen::Index npoints=static_cast<Eigen::Index>(angular.size());
    Eigen::Matrix<double,3,Eigen::Dynamic> points(3,npoints);
    Eigen::VectorXd weights(npoints), ones=Eigen::VectorXd::Ones(npoints);
    std::vector<double> edges{0.0,1e-5};
    for(int i=1;i<=16;++i)
        edges.push_back(1e-5*std::pow(4.0e6,double(i)/16.0));
    for(std::size_t a=0;a<centers.size();++a) {
        for(std::size_t panel=0;panel+1<edges.size();++panel) {
            const double lo=edges[panel],hi=edges[panel+1];
            for(int ir=0;ir<radial_order;++ir) {
                const double r=lo+(hi-lo)*radial.x[ir];
                const double wr=(hi-lo)*radial.w[ir]*r*r;
                for(Eigen::Index i=0;i<npoints;++i) {
                    const auto& v=angular[static_cast<std::size_t>(i)];
                    const Eigen::Vector3d p=centers[a].r+r*Eigen::Vector3d(v.x,v.y,v.z);
                    points.col(i)=p;
                    weights[i]=wr*v.w*becke_weights(p,centers)[static_cast<Eigen::Index>(a)];
                }
                accumulator.add(points,weights,ones);
            }
        }
    }
    return accumulator.matrix();
}
void validate_angular_rules() {
    for(int count : {194,302,434,590,770}) {
        const auto angular=lebedev_reference::grid(count);
        double zeroth=0.0, first=0.0, second=0.0, fourth=0.0;
        for(const auto& v:angular) {
            const double norm=v.x*v.x+v.y*v.y+v.z*v.z;
            if(std::abs(norm-1.0)>2e-13) throw std::runtime_error("Lebedev point off sphere");
            zeroth+=v.w;
            first+=v.w*v.x;
            second+=v.w*v.x*v.x;
            fourth+=v.w*v.x*v.x*v.y*v.y;
        }
        const double four_pi=4.0*std::numbers::pi;
        if(std::abs(zeroth-four_pi)>2e-12 || std::abs(first)>2e-12 ||
           std::abs(second-four_pi/3.0)>2e-12 ||
           std::abs(fourth-four_pi/15.0)>2e-12)
            throw std::runtime_error("Lebedev moment validation failed");
    }
}

}
int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        validate_angular_rules();
        for(int arg=1;arg<argc;++arg) {
            chem::kumo::FormattedCheckpoint fchk(argv[arg]);
            chem::kumo::BasisEvaluator ao{chem::kumo::BasisSet(fchk)};
            const auto centers=centers_of(ao.basis());
            const auto exact=ao.overlap_matrix();
            // Check partition of unity independently at representative off-nuclear points.
            double unity=0.0;
            for(const auto& c:centers) {
                const Eigen::Vector3d p=c.r+Eigen::Vector3d(0.31,-0.27,0.19);
                unity=std::max(unity,std::abs(becke_weights(p,centers).sum()-1.0));
            }
            if(unity>1e-13) throw std::runtime_error("Becke partition of unity failed");
            std::cout<<std::scientific<<std::setprecision(12)
                     <<argv[arg]<<" nAO="<<ao.basis_size()<<" centers="<<centers.size()
                     <<" partition_unity="<<unity;
            const auto atom_of_ao=ao_centers(ao.basis());
            // Fixed radial order, compare five Lebedev angular grids.
            constexpr int radial_order=20;
            for(const int angular_points : {194,302,434,590,770}) {
                const auto numerical=integrate_lebedev(ao,centers,radial_order,angular_points);
                const Eigen::MatrixXd difference=(numerical-exact).cwiseAbs();
                double intra=0.0, inter=0.0;
                for(Eigen::Index i=0;i<difference.rows();++i)
                    for(Eigen::Index j=0;j<difference.cols();++j) {
                        double& error=(atom_of_ao[static_cast<std::size_t>(i)]==
                                       atom_of_ao[static_cast<std::size_t>(j)])
                                      ? intra : inter;
                        error=std::max(error,difference(i,j));
                    }
                const double symmetry=(numerical-numerical.transpose()).cwiseAbs().maxCoeff();
                if(!std::isfinite(intra)||!std::isfinite(inter)||!std::isfinite(symmetry))
                    throw std::runtime_error("nonfinite molecular overlap");
                std::cout<<" radial="<<radial_order
                         <<" lebedev="<<angular_points
                         <<" intra="<<intra<<" inter="<<inter
                         <<" max_abs="<<std::max(intra,inter)
                         <<" symmetry="<<symmetry;
            }
            std::cout<<'\n';
        }
    } catch(const std::exception& e) {
        std::cerr<<"multicenter overlap diagnostic: "<<e.what()<<'\n';
        return 1;
    }
}
