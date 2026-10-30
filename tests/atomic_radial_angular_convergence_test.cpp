// Stage 03C.3: independent single-centre radial/angular quadrature validation.
// No AO evaluator or FCHK parsing: analytic Gaussian moments are the oracle.
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
struct Rule { std::vector<double> x, w; };
// Gauss-Legendre on [0,1], Golub-Welsch; independent of production AO code.
Rule legendre(int n) {
    Eigen::MatrixXd j = Eigen::MatrixXd::Zero(n,n);
    for (int k=1;k<n;++k) {
        const double b = double(k)/std::sqrt(4.0*k*k-1.0);
        j(k-1,k)=j(k,k-1)=b;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(j);
    if (eig.info()!=Eigen::Success) throw std::runtime_error("Legendre eigensolver failed");
    Rule r;
    for (int i=0;i<n;++i) {
        r.x.push_back((eig.eigenvalues()[i]+1.0)/2.0);
        r.w.push_back(std::pow(eig.eigenvectors()(0,i),2));
    }
    return r;
}
// Cartesian monomial moments of exp(-beta*r^2), integrated over all R^3.
// Monomials: 1, x^2, x^4, x^2*y^2, x^6.
std::array<double,5> exact(double beta) {
    const double z = std::pow(std::numbers::pi/beta,1.5);
    return {z,z/(2*beta),3*z/(4*beta*beta),z/(4*beta*beta),
            15*z/(8*beta*beta*beta)};
}
std::array<double,5> integrate(double beta, int nr, int na, int nphi) {
    const Rule radial=legendre(nr), polar=legendre(na);
    std::array<double,5> sums{};
    // r = scale*t/(1-t), t in (0,1). Jacobian r^2 dr/dt.
    const double scale=1.0/std::sqrt(beta);
    for (int ir=0;ir<nr;++ir) {
        const double t=radial.x[ir], r=scale*t/(1-t);
        const double rw=radial.w[ir]*r*r*scale/((1-t)*(1-t))*std::exp(-beta*r*r);
        for (int ia=0;ia<na;++ia) {
            const double mu=2*polar.x[ia]-1;
            const double sin_theta=std::sqrt(std::max(0.0,1-mu*mu));
            const double aw=2*polar.w[ia]*2*std::numbers::pi/nphi;
            for (int ip=0;ip<nphi;++ip) {
                const double phi=2*std::numbers::pi*(ip+0.5)/nphi;
                const double x=r*sin_theta*std::cos(phi);
                const double y=r*sin_theta*std::sin(phi);
                const double x2=x*x, y2=y*y, weight=rw*aw;
                sums[0]+=weight;
                sums[1]+=weight*x2;
                sums[2]+=weight*x2*x2;
                sums[3]+=weight*x2*y2;
                sums[4]+=weight*x2*x2*x2;
            }
        }
    }
    return sums;
}
}
int main() {
    try {
        // alpha ranges from diffuse to tight; beta=2*alpha for same-centre AO overlaps.
        for (double alpha : {0.075,0.5,10.0,100.0}) {
            const double beta=2*alpha;
            const auto ref=exact(beta);
            double previous=1e100;
            for (int nr : {12,24,48,64,80}) {
                const auto got=integrate(beta,nr,8,20);
                double worst=0;
                for (std::size_t k=0;k<got.size();++k)
                    worst=std::max(worst,std::abs(got[k]/ref[k]-1.0));
                std::cout<<std::scientific<<std::setprecision(9)
                         <<"alpha="<<alpha<<" radial="<<nr
                         <<" max_relative_moment_error="<<worst<<'\n';
                if (!std::isfinite(worst) || (nr==48 && worst>2e-9) || (nr==80 && worst>2e-11))
                    throw std::runtime_error("atomic radial-angular quadrature failed");
                if (nr>12 && worst>previous*1.05+1e-14)
                    throw std::runtime_error("quadrature error did not improve");
                previous=worst;
            }
        }
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
