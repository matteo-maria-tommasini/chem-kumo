// Exploratory 4-dimensional central-atom valence coupling diagnostic.
// No lone-pair labels are inferred from eigenvalues alone.
#include "chem/kumo/local_neighbor_w.hpp"
#include "chem/kumo/local_joint_diagonalization.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
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
    if(argc!=2){std::cerr<<"usage: kumo-tetrahedral-w-diagnostic ch4.fchk.gz\n";return 2;}
    try {
        std::cout<<std::setprecision(12);
        for(int arg=1;arg<argc;++arg) {

            const FormattedCheckpoint fchk(argv[arg]);
            const BasisSet basis(fchk);
            const Molecule mol=Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd p=mol.density_matrix(Spin::Total);
            const Eigen::MatrixXd s=BasisEvaluator(basis).overlap_matrix();
            const auto val=select_valence(basis);
            if(mol.atoms().size()!=5 || mol.atoms()[0].atomic_number!=6)
                throw std::runtime_error("expected CH4 with carbon as atom zero");
            Eigen::Matrix4d h;
            for(Eigen::Index i=0;i<4;++i) {
                if(mol.atoms()[static_cast<std::size_t>(i+1)].atomic_number!=1)
                    throw std::runtime_error("expected four hydrogens");
                const Eigen::Vector3d d=(mol.atoms()[static_cast<std::size_t>(i+1)].position-
                    mol.atoms()[0].position).normalized();
                h(0,i)=0.5;
                h.block<3,1>(1,i)=std::sqrt(3.0)*0.5*d;
            }
            const Eigen::Matrix4d gram=h.transpose()*h;
            const Eigen::Matrix4d ugeom=h*invsqrt(gram);
            std::cout<<"# geometric_hybrid_gram\n"<<gram<<"\n";
            std::cout<<"# geometric_hybrid_rotation_columns\n"<<ugeom<<"\n";
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
                per_neighbor.push_back(g);
            }
            std::cout<<"# W and G matrices in geometric hybrid basis; each neighbor indexed from 1\n";
            for(std::size_t b=1;b<mol.atoms().size();++b) {
                const auto ab=local_neighbor_w_block(basis,{p},0,b);
                const Eigen::MatrixXd xb=invsqrt(sub(s,ab.ao_b,ab.ao_b));
                std::vector<Eigen::Index> all;
                for(Eigen::Index j=0;j<ab.matrices[0].cols();++j) all.push_back(j);
                const Eigen::MatrixXd wab=sub(ab.matrices[0],local,all);
                const Eigen::MatrixXd t=(x*wab*xb).eval();
                const Eigen::MatrixXd transformed_t=(ugeom.transpose()*t).eval();
                const Eigen::Matrix4d transformed_g=ugeom.transpose()*per_neighbor[b-1]*ugeom;
                const double trace=transformed_g.trace();
                std::cout<<"# neighbor "<<b<<" normalized_W_hybrid_rows\n"<<transformed_t<<"\n";
                std::cout<<"# neighbor "<<b<<" G_hybrid\n"<<transformed_g<<"\n";
                std::cout<<"# neighbor "<<b<<" diagonal_fraction="
                         <<(trace>0.0?transformed_g(b-1,b-1)/trace:0.0)<<"\n";
            }
            std::cerr<<"system="<<argv[arg]<<" neighbors=4 status=OK\n";
        }
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
