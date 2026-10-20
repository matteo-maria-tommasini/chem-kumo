// Stage 03C.8b.1: independent AO-evaluator quadrature for contracted quartets.
#include "chem/kumo/local_shell_integrals.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
namespace {
struct Rule { Eigen::VectorXd nodes, weights; };
Rule legendre(int n) {
    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(n,n);
    for (int k=1;k<n;++k) {
        const double b=double(k)/std::sqrt(4.0*k*k-1.0);
        jac(k-1,k)=jac(k,k-1)=b;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(jac);
    if (eig.info()!=Eigen::Success) throw std::runtime_error("Legendre diagonalization failed");
    return {(eig.eigenvalues().array()+1.0).matrix()/2.0,
            eig.eigenvectors().row(0).array().square().transpose().matrix()};
}
struct Case { std::array<std::size_t,4> shell; const char* name; };
// A single global spherical integration, with no Becke weights. This tests
// the contracted shell kernel against the *independent* AO evaluator.
double quadrature(const BasisEvaluator& ao, const std::array<Eigen::Index,4>& idx,
                  const Eigen::Vector3d& origin, int radial_order) {
    constexpr int nmu=24, nphi=48;
    const auto radial=legendre(radial_order), polar=legendre(nmu);
    std::vector<double> edges{0.0,1e-5};
    for (int i=1;i<=16;++i)
        edges.push_back(1e-5*std::pow(4.0e6,double(i)/16.0));
    double result=0.0;
    for (std::size_t panel=0;panel+1<edges.size();++panel) {
        const double lo=edges[panel], hi=edges[panel+1];
        for (int ir=0;ir<radial_order;++ir) {
            const double r=lo+(hi-lo)*radial.nodes[ir];
            const double wr=(hi-lo)*radial.weights[ir]*r*r;
            for (int im=0;im<nmu;++im) {
                const double mu=2.0*polar.nodes[im]-1.0;
                const double st=std::sqrt(std::max(0.0,1.0-mu*mu));
                for (int ip=0;ip<nphi;++ip) {
                    const double phi=2.0*std::numbers::pi*(ip+0.5)/nphi;
                    const Eigen::Vector3d pos=origin+r*Eigen::Vector3d(st*std::cos(phi),st*std::sin(phi),mu);
                    const auto values=ao.evaluate(pos);
                    result += wr*2.0*polar.weights[im]*(2.0*std::numbers::pi/nphi)
                        *values[idx[0]]*values[idx[1]]*values[idx[2]]*values[idx[3]];
                }
            }
        }
    }
    return result;
}
}
int main(int argc,char** argv) {
    if (argc<2) return 2;
    try {
        for (int arg=1;arg<argc;++arg) {
            const FormattedCheckpoint fchk(argv[arg]);
            const BasisSet basis(fchk);
            const BasisEvaluator ao(basis);
            const auto components=local_shell_components(basis);
            const auto& shells=basis.shells();
            std::vector<Eigen::Index> offsets;
            std::vector<int> angular;
            std::vector<std::size_t> centers;
            Eigen::Index offset=0;
            for (const auto& sh:shells) {
                offsets.push_back(offset);
                if (sh.is_sp) {
                    angular.push_back(0); centers.push_back(sh.center);
                    offsets.push_back(offset+1);
                    angular.push_back(1); centers.push_back(sh.center);
                } else {
                    angular.push_back(sh.angular_momentum); centers.push_back(sh.center);
                }
                offset+=static_cast<Eigen::Index>(sh.function_count());
            }
            if (offset!=static_cast<Eigen::Index>(ao.basis_size()) || offsets.size()!=components.size())
                throw std::runtime_error("shell/AO offsets mismatch");
            const auto find_l=[&](int l) -> std::size_t {
                for (std::size_t i=0;i<angular.size();++i)
                    if (angular[i]==l) return i;
                return angular.size();
            };
            const auto s=find_l(0), p=find_l(1), d=find_l(2);
            if (s==components.size() || p==components.size())
                throw std::runtime_error("test requires s and p shells");
            std::vector<Case> cases{{{s,s,s,s},"ssss"},{{p,s,p,s},"psps"}};
            if (d<components.size()) cases.push_back({{d,s,d,s},"dsds"});
            // Exercise a genuinely multicentric quartet where available.
            for (std::size_t j=0;j<components.size();++j)
                if (centers[j]!=centers[s]) {cases.push_back({{s,j,s,j},"two-center"});break;}
            double maximum=0.0;
            for (const auto& c:cases) {
                const auto block=local_shell_quartet({&components[c.shell[0]],&components[c.shell[1]],
                                                       &components[c.shell[2]],&components[c.shell[3]]});
                const std::array<Eigen::Index,4> idx{offsets[c.shell[0]],offsets[c.shell[1]],
                                                      offsets[c.shell[2]],offsets[c.shell[3]]};
                Eigen::Vector3d center=Eigen::Vector3d::Zero();
                for (const auto i:c.shell) {
                    // All split SP parts share the same original center.
                    const auto original_center=centers[i];
                    const auto it=std::find_if(shells.begin(),shells.end(),
                        [&](const Shell& sh){return sh.center==original_center;});
                    if (it==shells.end()) throw std::runtime_error("missing center");
                    center+=it->origin/4.0;
                }
                const double exact=block.at(0,0,0,0);
                const double q12=quadrature(ao,idx,center,12);
                const double q16=quadrature(ao,idx,center,16);
                const double err=std::abs(exact-q16);
                maximum=std::max(maximum,err);
                std::cout<<std::scientific<<std::setprecision(12)<<argv[arg]
                         <<" case="<<c.name<<" analytic="<<exact
                         <<" q12="<<q12<<" q16="<<q16<<" abs_error="<<err<<'\n';
                if (!std::isfinite(err) || err>2e-7*std::max(1.0,std::abs(exact)))
                    throw std::runtime_error("contracted quartet numerical reference mismatch");
            }
            std::cout<<argv[arg]<<" checked="<<cases.size()<<" max_abs="<<maximum<<'\n';
        }
    } catch(const std::exception& e) {
        std::cerr<<"local shell quartet reference: "<<e.what()<<'\n'; return 1;
    }
}
