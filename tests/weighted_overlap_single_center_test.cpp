// Stage 03C.5a: end-to-end, all atomic-center AO overlaps by radial/angular quadrature.
// The integration rule is deliberately test-only; no production quadrature API.
#include "chem/kumo/weighted_overlap.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
struct Rule { Eigen::VectorXd x, w; };
Rule legendre(int n) {
    Eigen::MatrixXd j = Eigen::MatrixXd::Zero(n,n);
    for (int k=1;k<n;++k) {
        const double b=double(k)/std::sqrt(4.0*k*k-1.0);
        j(k-1,k)=j(k,k-1)=b;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(j);
    if(eig.info()!=Eigen::Success) throw std::runtime_error("Gauss-Legendre eigensolver failed");
    return {(eig.eigenvalues().array()+1.0).matrix()/2.0,
            eig.eigenvectors().row(0).array().square().transpose().matrix()};
}

// Logarithmically spaced radial panels cover multiple primitive length scales.
// Every panel uses independent Gauss-Legendre nodes (no adaptive feedback).
Eigen::MatrixXd integrate(const chem::kumo::BasisEvaluator& ao,
                          const Eigen::Vector3d& origin,
                          int radial_order) {
    constexpr int nmu=10, nphi=24;
    const auto radial=legendre(radial_order), polar=legendre(nmu);
    chem::kumo::WeightedOverlapAccumulator accumulator(ao);
    Eigen::Matrix<double,3,Eigen::Dynamic> points(3,nmu*nphi);
    Eigen::VectorXd weights(nmu*nphi), ones=Eigen::VectorXd::Ones(nmu*nphi);
    // 16 logarithmic panels from 1e-5 to 40 bohr, plus [0,1e-5].
    std::vector<double> edges{0.0,1e-5};
    for(int i=1;i<=16;++i)
        edges.push_back(1e-5*std::pow(4.0e6, double(i)/16.0));
    for(std::size_t panel=0;panel+1<edges.size();++panel) {
        const double lo=edges[panel], hi=edges[panel+1];
        for(int ir=0;ir<radial_order;++ir) {
            const double r=lo+(hi-lo)*radial.x[ir];
            const double wr=(hi-lo)*radial.w[ir]*r*r;
            for(int im=0;im<nmu;++im) {
                const double mu=2*polar.x[im]-1.0;
                const double sint=std::sqrt(std::max(0.0,1.0-mu*mu));
                for(int ip=0;ip<nphi;++ip) {
                    const double phi=2*std::numbers::pi*(ip+0.5)/nphi;
                    const Eigen::Index idx=im*nphi+ip;
                    points.col(idx)=origin+r*Eigen::Vector3d(
                        sint*std::cos(phi),sint*std::sin(phi),mu);
                    weights[idx]=wr*2*polar.w[im]*(2*std::numbers::pi/nphi);
                }
            }
            accumulator.add(points,weights,ones);
        }
    }
    return accumulator.matrix();
}

struct AtomicCenter {
    std::size_t center;
    Eigen::Vector3d origin;
    std::vector<Eigen::Index> indices;
};

std::vector<AtomicCenter> atomic_centers(const chem::kumo::BasisSet& basis) {
    const auto& shells=basis.shells();
    if(shells.empty()) throw std::runtime_error("empty basis");
    std::vector<AtomicCenter> centers;
    Eigen::Index offset=0;
    for(const auto& shell:shells) {
        auto it=std::find_if(centers.begin(),centers.end(),
            [&](const AtomicCenter& c) { return c.center==shell.center; });
        if(it==centers.end()) {
            centers.push_back({static_cast<std::size_t>(shell.center),
                               shell.origin,{}});
            it=std::prev(centers.end());
        } else if((it->origin-shell.origin).norm()>1e-10) {
            throw std::runtime_error("inconsistent origins for an atomic center");
        }
        for(std::size_t j=0;j<shell.function_count();++j)
            it->indices.push_back(offset+static_cast<Eigen::Index>(j));
        offset+=static_cast<Eigen::Index>(shell.function_count());
    }
    if(offset!=static_cast<Eigen::Index>(basis.function_count()))
        throw std::runtime_error("AO shell offsets inconsistent with basis size");
    return centers;
}

double block_error(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b,
                   const std::vector<Eigen::Index>& indices) {
    double error=0.0;
    for(auto i:indices) for(auto j:indices)
        error=std::max(error,std::abs(a(i,j)-b(i,j)));
    return error;
}
}

int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        for(int arg=1;arg<argc;++arg) {
            chem::kumo::FormattedCheckpoint fchk(argv[arg]);
            chem::kumo::BasisEvaluator ao{chem::kumo::BasisSet(fchk)};
            const auto centers=atomic_centers(ao.basis());
            const auto exact=ao.overlap_matrix();
            for(const auto& center:centers) {
                std::cout<<std::scientific<<std::setprecision(12)
                         <<argv[arg]<<" nAO="<<ao.basis().function_count()
                         <<" center="<<center.center
                         <<" center_AOs="<<center.indices.size();
                double previous=0.0;
                for(int n : {8,12,16}) {
                    const auto numeric=integrate(ao,center.origin,n);
                    const double error=block_error(numeric,exact,center.indices);
                    if(!std::isfinite(error))
                        throw std::runtime_error("nonfinite overlap error");
                    std::cout<<" radial_per_panel="<<n<<" error="<<error;
                    if(n==16 && error>3e-7)
                        throw std::runtime_error("atomic AO overlap did not converge");
                    if(n==16 && error>previous*1.1+1e-11)
                        throw std::runtime_error("atomic AO overlap worsened");
                    previous=error;
                }
                std::cout<<'\n';
            }
        }
    } catch(const std::exception& e) {
        std::cerr<<"single-center AO quadrature: "<<e.what()<<'\n';
        return 1;
    }
}
