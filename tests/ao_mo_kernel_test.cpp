#include "chem/kumo/density.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace chem::kumo::detail {
struct AoMoTestAccess {
    static auto run(const Molecule& m,
        const Eigen::Ref<const Eigen::Matrix<double,3,Eigen::Dynamic>>& points,
        const std::vector<Spin>& components, std::size_t block, int order) {
        return m.density_derivatives_mo(points,components,block,order);
    }
    static auto reconstructed_matrices(const Molecule& m) {
        const auto& cb = m.unrestricted_ ? m.beta_mos_ : m.alpha_mos_;
        const Eigen::MatrixXd pa = m.alpha_mos_.leftCols(m.nalpha_) *
                                   m.alpha_mos_.leftCols(m.nalpha_).transpose();
        const Eigen::MatrixXd pb = cb.leftCols(m.nbeta_) * cb.leftCols(m.nbeta_).transpose();
        std::map<Spin,Eigen::MatrixXd> result;
        result.emplace(Spin::Alpha,pa);
        result.emplace(Spin::Beta,pb);
        result.emplace(Spin::Total,pa+pb);
        result.emplace(Spin::Difference,pa-pb);
        return result;
    }
    static DensityDerivatives ao_p_reference(const Molecule& m,
                                             const Eigen::Vector3d& point,
                                             const Eigen::MatrixXd& p) {
        const auto ao=m.evaluator_.evaluate_derivatives(point,2);
        const auto y=p*ao.values;
        const auto pg=p*ao.gradient;
        constexpr int index[3][3]={{0,3,4},{3,1,5},{4,5,2}};
        DensityDerivatives d;
        d.value=ao.values.dot(y);
        for(int a=0;a<3;++a) {
            d.gradient[a]=2.0*ao.gradient.col(a).dot(y);
            d.tau+=0.5*ao.gradient.col(a).dot(pg.col(a));
            for(int b=a;b<3;++b)
                d.hessian(a,b)=d.hessian(b,a)=2.0*(
                    ao.hessian.col(index[a][b]).dot(y)+
                    ao.gradient.col(a).dot(pg.col(b)));
        }
        d.laplacian=d.hessian.trace();
        return d;
    }
};
}
using namespace chem::kumo;
namespace {
constexpr double algorithm_tol=2e-10;
void check(double actual,double reference,const std::string& label,
           const char* file,Spin spin,std::size_t block,std::size_t point) {
    const double scale=std::max({1.0,std::abs(actual),std::abs(reference)});
    const double error=std::abs(actual-reference);
    if(!std::isfinite(actual) || !std::isfinite(reference) || error>algorithm_tol*scale) {
        std::ostringstream msg;
        msg<<std::scientific<<std::setprecision(16)
           <<"AO-MO algorithm mismatch: "<<label
           <<" file="<<file<<" spin="<<static_cast<int>(spin)
           <<" block="<<block<<" point="<<point
           <<"\n  AO-MO="<<actual<<"\n  AO-P(P_MO)="<<reference
           <<"\n  abs_error="<<error<<"\n  scaled_error="<<error/scale
           <<"\n  tolerance="<<algorithm_tol*scale;
        throw std::runtime_error(msg.str());
    }
}
void compare(const DensityDerivatives& x,const DensityDerivatives& y,
             const char* file,Spin spin,std::size_t block,std::size_t point,int order) {
    const auto test=[&](double a,double b,const std::string& name) {
        check(a,b,name,file,spin,block,point);
    };
    test(x.value,y.value,"density");
    for(int a=0;a<3;++a) test(x.gradient[a],y.gradient[a],std::string("gradient_")+"xyz"[a]);
    test(x.tau,y.tau,"tau");
    if(order==1) {
        test(x.laplacian,0.0,"first-order laplacian");
        return;
    }
    test(x.laplacian,y.laplacian,"laplacian");
    constexpr int pairs[6][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
    constexpr const char* names[6]={"xx","yy","zz","xy","xz","yz"};
    for(int k=0;k<6;++k)
        test(x.hessian(pairs[k][0],pairs[k][1]),
             y.hessian(pairs[k][0],pairs[k][1]),std::string("hessian_")+names[k]);
    test(x.hessian.trace(),x.laplacian,"MO Hessian trace");
    test(y.hessian.trace(),y.laplacian,"AO-P Hessian trace");
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_mo(x.hessian);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_ao(y.hessian);
    if(eig_mo.info()!=Eigen::Success || eig_ao.info()!=Eigen::Success)
        throw std::runtime_error("Hessian eigensolver failed");
    const double bound=(x.hessian-y.hessian).norm(); // Weyl: |delta lambda_i| <= ||delta H||_2 <= ||delta H||_F
    for(int k=0;k<3;++k) {
        const double diff=std::abs(eig_mo.eigenvalues()[k]-eig_ao.eigenvalues()[k]);
        test(eig_mo.eigenvalues()[k],eig_ao.eigenvalues()[k],
             "hessian_eigenvalue_"+std::to_string(k));
        if(diff>bound+1e-10*std::max({1.0,eig_mo.eigenvalues().cwiseAbs().maxCoeff(),
                                      eig_ao.eigenvalues().cwiseAbs().maxCoeff()}))
            throw std::runtime_error("Weyl bound violated");
    }
    const double lambda2_mo=eig_mo.eigenvalues()[1];
    const double lambda2_ao=eig_ao.eigenvalues()[1];
    // Only compare signs when both are well separated from the error bound.
    if(std::min(std::abs(lambda2_mo),std::abs(lambda2_ao))>bound+1e-10 &&
       std::signbit(lambda2_mo)!=std::signbit(lambda2_ao))
        throw std::runtime_error("Stable lambda2 sign mismatch");
}
}
int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        const std::vector<Spin> spins{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference};
        Eigen::Matrix<double,3,Eigen::Dynamic> points(3,19);
        for(int j=0;j<19;++j)
            points.col(j)=Eigen::Vector3d(0.17*j-1.2,0.09*j-0.8,0.13*j+0.2);
        for(int arg=1;arg<argc;++arg) {
            const auto m=Molecule::from_fchk(argv[arg]);
            if(!m.electronic_structure_info().ao_mo_eligible())
                throw std::runtime_error("test molecule not AO-MO eligible");
            const auto matrices=detail::AoMoTestAccess::reconstructed_matrices(m);
            for(int order:{1,2})
                for(std::size_t block:{std::size_t{1},std::size_t{7},std::size_t{32}}) {
                    const auto mo=detail::AoMoTestAccess::run(m,points,spins,block,order);
                    for(Spin spin:spins)
                        for(std::size_t j=0;j<static_cast<std::size_t>(points.cols());++j) {
                            const auto ref=detail::AoMoTestAccess::ao_p_reference(
                                m,points.col(static_cast<Eigen::Index>(j)),matrices.at(spin));
                            compare(mo.at(spin).at(j),ref,argv[arg],spin,block,j,order);
                        }
                }
        }
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
