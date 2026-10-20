#include "chem/kumo/local_gaussian_integrals.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <numbers>
#include <string>

using chem::kumo::LocalGaussianPrimitive;
using chem::kumo::local_four_gaussian_integral;

// Independent tensor-product Gauss-Hermite integration. The Gaussian product
// center and exponent are used only to choose the integration coordinates;
// polynomial factors are evaluated directly at each quadrature node.
static double gauss_hermite_reference(const std::array<LocalGaussianPrimitive,4>& g) {
    constexpr int n = 16; // exact through polynomial degree 31 in each axis
    Eigen::MatrixXd jacobi = Eigen::MatrixXd::Zero(n,n);
    for (int i=0;i<n-1;++i) {
        const double v=std::sqrt((i+1)/2.0);
        jacobi(i,i+1)=v;
        jacobi(i+1,i)=v;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(jacobi);
    if (solver.info()!=Eigen::Success) throw std::runtime_error("Hermite eigensolver failed");
    std::array<double,n> nodes{}, weights{};
    for(int i=0;i<n;++i) {
        nodes[i]=solver.eigenvalues()[i];
        weights[i]=std::sqrt(std::numbers::pi)*std::pow(solver.eigenvectors()(0,i),2);
    }
    double p=0, scale=1;
    std::array<double,3> center{};
    for(const auto& a:g) {
        p+=a.exponent;
        scale*=a.scale;
        for(int d=0;d<3;++d) center[d]+=a.exponent*a.center[d];
    }
    for(auto& x:center) x/=p;
    double separation=0;
    for(std::size_t i=0;i<4;++i) for(std::size_t j=i+1;j<4;++j) {
        double r2=0;
        for(int d=0;d<3;++d) r2+=std::pow(g[i].center[d]-g[j].center[d],2);
        separation+=g[i].exponent*g[j].exponent/p*r2;
    }
    const double inv_sqrt_p=1/std::sqrt(p);
    double sum=0;
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) for(int k=0;k<n;++k) {
        const std::array<double,3> r{
            center[0]+nodes[i]*inv_sqrt_p,
            center[1]+nodes[j]*inv_sqrt_p,
            center[2]+nodes[k]*inv_sqrt_p};
        double polynomial=1;
        for(const auto& a:g) for(int d=0;d<3;++d)
            polynomial*=std::pow(r[d]-a.center[d],a.powers[d]);
        sum+=weights[i]*weights[j]*weights[k]*polynomial;
    }
    return scale*std::exp(-separation)*sum*inv_sqrt_p*inv_sqrt_p*inv_sqrt_p;
}

static void compare(const std::array<LocalGaussianPrimitive,4>& g, const char* name) {
    const double analytic=local_four_gaussian_integral(g);
    const double reference=gauss_hermite_reference(g);
    const double error=std::abs(analytic-reference);
    const double tolerance=2e-11*std::max({1.0,std::abs(analytic),std::abs(reference)});
    std::cout<<name<<" analytic="<<std::setprecision(16)<<analytic
             <<" reference="<<reference<<" error="<<error<<'\n';
    if(!std::isfinite(analytic)||error>tolerance)
        throw std::runtime_error(std::string("four-Gaussian mismatch: ")+name);
}

int main() {
    try {
        LocalGaussianPrimitive a{{0,0,0},{0,0,0},0.7,1};
        LocalGaussianPrimitive b{{0.3,-0.2,0.5},{0,0,0},1.2,0.8};
        LocalGaussianPrimitive c{{-0.6,0.4,0.2},{0,0,0},0.35,-0.7};
        LocalGaussianPrimitive d{{0.1,0.8,-0.3},{0,0,0},2.4,1.3};
        compare({a,b,c,d},"ssss");
        a.powers={1,0,0}; compare({a,b,c,d},"psss");
        b.powers={0,1,0}; compare({a,b,c,d},"ppss");
        c.powers={2,0,0}; compare({a,b,c,d},"ppds");
        d.powers={1,1,1}; compare({a,b,c,d},"ppdf");
        a.powers={3,0,0}; b.powers={0,2,0};
        c.powers={1,0,2}; d.powers={0,1,2};
        compare({a,b,c,d},"high-angular-momentum");
        a.center={0,0,0}; b.center={0,0,0};
        c.center={0,0,0}; d.center={0,0,0};
        compare({a,b,c,d},"same-center");
        bool rejected=false;
        a.exponent=-1;
        try { (void)local_four_gaussian_integral({a,b,c,d}); }
        catch(const std::invalid_argument&) { rejected=true; }
        if(!rejected) throw std::runtime_error("negative exponent accepted");
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
