// Exploratory 4-dimensional central-atom valence coupling diagnostic.
// No lone-pair labels are inferred from eigenvalues alone.
#include "chem/kumo/local_neighbor_w.hpp"
#include "chem/kumo/local_joint_diagonalization.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;

static Eigen::MatrixXd sub(const Eigen::MatrixXd& m,
                           const std::vector<Eigen::Index>& a,
                           const std::vector<Eigen::Index>& b) {
    Eigen::MatrixXd out(a.size(),b.size());
    for(Eigen::Index i=0;i<out.rows();++i)
        for(Eigen::Index j=0;j<out.cols();++j)
            out(i,j)=m(a[static_cast<std::size_t>(i)],b[static_cast<std::size_t>(j)]);
    return out;
}

static Eigen::MatrixXd invsqrt(const Eigen::MatrixXd& s) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
    if(es.info()!=Eigen::Success || es.eigenvalues().minCoeff()<=1e-10*es.eigenvalues().maxCoeff())
        throw std::runtime_error("local overlap not positive definite at requested tolerance");
    Eigen::MatrixXd v=es.eigenvectors();
    Eigen::VectorXd d=es.eigenvalues().array().rsqrt();
    return (v*d.asDiagonal()*v.transpose()).eval();
}

// Select the second tightest contracted s shell and the tightest contracted p shell.
// This is an EXPLICIT, basis-dependent exploratory convention, not a canonical
// definition of the atomic valence subspace.
static std::vector<Eigen::Index> select_valence(const BasisSet& basis) {
    struct Choice { double exponent; std::vector<Eigen::Index> ao; };
    std::vector<Choice> ss,ps;
    Eigen::Index offset=0;
    for(const auto& shell:basis.shells()) {
        const auto count=static_cast<Eigen::Index>(shell.function_count());
        if(shell.center==0) {
            double e=0.0;
            for(double x:shell.exponents)e=std::max(e,x);
            if(shell.is_sp) {
                ss.push_back({e,{offset}});
                ps.push_back({e,{offset+1,offset+2,offset+3}});
            } else if(shell.angular_momentum==0 && count==1) {
                ss.push_back({e,{offset}});
            } else if(shell.angular_momentum==1 && count==3) {
                ps.push_back({e,{offset,offset+1,offset+2}});
            }
        }
        offset+=count;
    }
    auto descending=[](const Choice& a,const Choice& b){return a.exponent>b.exponent;};
    std::sort(ss.begin(),ss.end(),descending);
    std::sort(ps.begin(),ps.end(),descending);
    if(ss.size()<2 || ps.empty())throw std::runtime_error("no unambiguous s/p valence shell candidates");
    auto indices=ss[1].ao;
    indices.insert(indices.end(),ps[0].ao.begin(),ps[0].ao.end());
    return indices;
}

int main(int argc,char** argv) {
    if(argc!=4){std::cerr<<"usage: diagnostic h2o.fchk.gz nh3.fchk.gz ch4.fchk.gz\n";return 2;}
    try {
        std::cout<<std::setprecision(12);
        // Mixed diagnostic records: AJD rows include one coupling per neighbor; legacy rows follow.
        std::cout<<"# ajd,system,mode,total_coupling,per_neighbor_couplings...,ao_s,ao_px,ao_py,ao_pz\n";
        for(int arg=1;arg<argc;++arg) {
            const FormattedCheckpoint fchk(argv[arg]);
            const BasisSet basis(fchk);
            const Molecule mol=Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd p=mol.density_matrix(Spin::Total);
            const Eigen::MatrixXd s=BasisEvaluator(basis).overlap_matrix();
            const auto val=select_valence(basis);
            const Eigen::MatrixXd saa=sub(s,val,val);
            const Eigen::MatrixXd x=invsqrt(saa);
            const auto aa=local_atomic_w(basis,{p},0);
            // Map selected global AO indices into the atomic W block.
            std::vector<Eigen::Index> local;
            for(auto gi:val){
                auto it=std::find(aa.ao_indices.begin(),aa.ao_indices.end(),gi);
                if(it==aa.ao_indices.end())throw std::runtime_error("valence AO absent from local block");
                local.push_back(static_cast<Eigen::Index>(it-aa.ao_indices.begin()));
            }
            const Eigen::MatrixXd w=sub(aa.matrices.at(0),local,local);
            const Eigen::MatrixXd pd=sub(p,val,val);
            const Eigen::MatrixXd wlocal=(x*w*x).eval();
            Eigen::MatrixXd m=Eigen::MatrixXd::Zero(4,4);
            std::vector<Eigen::MatrixXd> per_neighbor;
            for(std::size_t b=1;b<mol.atoms().size();++b) {
                const auto ab=local_neighbor_w_block(basis,{p},0,b);
                const Eigen::MatrixXd sab=sub(s,ab.ao_b,ab.ao_b);
                const Eigen::MatrixXd xb=invsqrt(sab);
                // AB target rows follow aa.ao_indices order.
                const Eigen::MatrixXd wab=sub(ab.matrices.at(0),local,
                    [&](){std::vector<Eigen::Index> ids;for(Eigen::Index j=0;j<ab.matrices[0].cols();++j)ids.push_back(j);return ids;}());
                const Eigen::MatrixXd t=(x*wab*xb).eval();
                const Eigen::MatrixXd g=(t*t.transpose()).eval();
                m+=g;
                per_neighbor.push_back(g);
            }
            const auto ajd=local_joint_diagonalize(per_neighbor);
            std::cerr<<"AJD system="<<argv[arg]<<" initial_offdiag="<<ajd.initial_offdiag
                     <<" final_offdiag="<<ajd.final_offdiag<<" sweeps="<<ajd.sweeps<<"\n";
            for(Eigen::Index i=0;i<4;++i) {
                const Eigen::VectorXd u=ajd.rotation.col(i);
                const Eigen::VectorXd c=x*u;
                std::cout<<"ajd,"<<argv[arg]<<","<<i<<","<<(u.transpose()*m*u)(0,0);
                for(const auto& g:ajd.transformed)std::cout<<","<<g(i,i);
                for(Eigen::Index j=0;j<4;++j)std::cout<<","<<c[j];
                std::cout<<"\n";
            }
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m);
            if(es.info()!=Eigen::Success)throw std::runtime_error("coupling eigensolver failed");
            for(Eigen::Index i=0;i<4;++i) {
                const Eigen::VectorXd v=es.eigenvectors().col(i);
                const Eigen::VectorXd c=x*v;
                const double dw=(c.transpose()*w*c)(0,0);
                const double dp=(c.transpose()*pd*c)(0,0); // NOT an electron occupation
                std::cout<<argv[arg]<<",4,"<<i<<","<<es.eigenvalues()[i]<<","<<dw<<","<<dp;
                for(Eigen::Index j=0;j<4;++j)std::cout<<","<<c[j];
                std::cout<<"\n";
            }
            std::cerr<<"system="<<argv[arg]<<" neighbors="<<(mol.atoms().size()-1)
                     <<" selected_global_AO="<<val[0]<<","<<val[1]<<","<<val[2]<<","<<val[3]
                     <<" status=OK\n";
        }
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
