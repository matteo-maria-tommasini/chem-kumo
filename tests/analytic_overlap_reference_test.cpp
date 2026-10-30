#include "chem/kumo/basis_evaluator.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
namespace {
struct Block {
    const Shell* shell;
    int l;
    std::vector<double> coefficients;
    std::vector<CartesianPower> powers;
    Eigen::MatrixXd transform;
    Eigen::Index offset;
};

// Independent Gaussian-Hermite quadrature: integrates exp(-x*x) times
// polynomials exactly through degree 2*n-1 (up to roundoff).
struct HermiteRule {
    Eigen::VectorXd nodes, weights;
    explicit HermiteRule(int n) {
        Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(n,n);
        for (int i=1;i<n;++i) jac(i-1,i)=jac(i,i-1)=std::sqrt(i/2.0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(jac);
        if (solver.info()!=Eigen::Success) throw std::runtime_error("Hermite rule failed");
        nodes=solver.eigenvalues();
        weights=(solver.eigenvectors().row(0).array().square()*std::sqrt(std::acos(-1.0))).matrix().transpose();
    }
};

double power(double x,int n) { double y=1; for(int i=0;i<n;++i)y*=x; return y; }
double odd_double_factorial(int n) { double v=1; for(int i=n;i>1;i-=2)v*=i; return v; }
double norm(double alpha,CartesianPower p) {
    const int l=p.x+p.y+p.z;
    const double pi=std::acos(-1.0);
    return std::pow(2*alpha/pi,0.75)*std::sqrt(power(4*alpha,l)/
        (odd_double_factorial(2*p.x-1)*odd_double_factorial(2*p.y-1)*odd_double_factorial(2*p.z-1)));
}

double axis_quadrature(int a,int b,double alpha,double beta,double A,double B,const HermiteRule& rule) {
    const double p=alpha+beta;
    const double center=(alpha*A+beta*B)/p;
    double sum=0;
    for (Eigen::Index k=0;k<rule.nodes.size();++k) {
        const double x=center+rule.nodes[k]/std::sqrt(p);
        sum+=rule.weights[k]*power(x-A,a)*power(x-B,b);
    }
    return std::exp(-alpha*beta/p*(A-B)*(A-B))*sum/std::sqrt(p);
}

std::vector<Block> blocks(const BasisSet& basis) {
    std::vector<Block> result;
    Eigen::Index offset=0;
    for(const auto& s:basis.shells()) {
        const auto add=[&](int l,const std::vector<double>& coeff,ShellRepresentation repr) {
            Block b{&s,l,coeff,gaussian_cartesian_order(l),{},offset};
            if(repr==ShellRepresentation::Pure && l>=2) b.transform=gaussian_cartesian_to_pure(l);
            offset+=b.transform.size()?b.transform.rows():static_cast<Eigen::Index>(b.powers.size());
            result.push_back(std::move(b));
        };
        if(s.is_sp) {add(0,s.coefficients,ShellRepresentation::Cartesian);add(1,s.p_coefficients,ShellRepresentation::Cartesian);}
        else add(s.angular_momentum,s.coefficients,s.representation);
    }
    if(offset!=static_cast<Eigen::Index>(basis.function_count())) throw std::runtime_error("AO offsets mismatch");
    return result;
}

Eigen::MatrixXd reference_block(const Block& a,const Block& b,const HermiteRule& rule) {
    Eigen::MatrixXd cart(a.powers.size(),b.powers.size());
    for(Eigen::Index i=0;i<cart.rows();++i) for(Eigen::Index j=0;j<cart.cols();++j) {
        const auto x=a.powers[static_cast<std::size_t>(i)];
        const auto y=b.powers[static_cast<std::size_t>(j)];
        double sum=0;
        for(std::size_t p=0;p<a.shell->exponents.size();++p)
            for(std::size_t q=0;q<b.shell->exponents.size();++q) {
                const double alpha=a.shell->exponents[p], beta=b.shell->exponents[q];
                double v=a.coefficients[p]*b.coefficients[q]*norm(alpha,x)*norm(beta,y);
                v*=axis_quadrature(x.x,y.x,alpha,beta,a.shell->origin.x(),b.shell->origin.x(),rule);
                v*=axis_quadrature(x.y,y.y,alpha,beta,a.shell->origin.y(),b.shell->origin.y(),rule);
                v*=axis_quadrature(x.z,y.z,alpha,beta,a.shell->origin.z(),b.shell->origin.z(),rule);
                sum+=v;
            }
        cart(i,j)=sum;
    }
    if(a.transform.size()) cart=(a.transform*cart).eval();
    if(b.transform.size()) cart=(cart*b.transform.transpose()).eval();
    return cart;
}
}
int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        const HermiteRule rule(12);
        for(int file=1;file<argc;++file) {
            FormattedCheckpoint fchk(argv[file]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            const auto S=evaluator.overlap_matrix();
            const auto bs=blocks(evaluator.basis());
            // Select diverse angular momenta, SP components and separated centers.
            // Bound runtime for large HBC fixtures; every selected block is checked fully.
            std::vector<std::size_t> selected;
            for(std::size_t i=0;i<bs.size();++i) {
                bool new_type=true;
                for(auto j:selected) if(bs[j].l==bs[i].l && bs[j].shell->center==bs[i].shell->center) {new_type=false;break;}
                if(new_type && selected.size()<12) selected.push_back(i);
            }
            double worst=0;
            std::size_t compared=0;
            for(auto ia:selected) for(auto ib:selected) {
                const auto& a=bs[ia]; const auto& b=bs[ib];
                const auto ref=reference_block(a,b,rule);
                const auto got=S.block(a.offset,b.offset,ref.rows(),ref.cols());
                const double error=(ref-got).cwiseAbs().maxCoeff();
                worst=std::max(worst,error);
                compared+=static_cast<std::size_t>(ref.size());
                if(!std::isfinite(error) || error>2e-11)
                    throw std::runtime_error("Gaussian-Hermite overlap reference mismatch");
            }
            std::cout<<argv[file]<<" reference_pairs="<<compared<<" max_abs="
                     <<std::scientific<<std::setprecision(12)<<worst<<'\n';
        }
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
