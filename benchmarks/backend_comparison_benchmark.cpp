#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Matrix = Eigen::MatrixXd;
using chem::kumo::FormattedCheckpoint;

Matrix unpack(const std::vector<double>& packed, Eigen::Index n) {
    if (packed.size() != static_cast<std::size_t>(n*(n+1)/2))
        throw std::runtime_error("invalid density matrix length");
    Matrix p(n,n);
    std::size_t k=0;
    for (Eigen::Index i=0;i<n;++i)
        for (Eigen::Index j=0;j<=i;++j) { p(i,j)=packed[k]; p(j,i)=packed[k]; ++k; }
    return p;
}
Matrix coefficients(const FormattedCheckpoint& f, const std::string& key, Eigen::Index n,
                    Eigen::Index occupied) {
    if (!f.contains(key)) throw std::runtime_error("missing " + key);
    const auto& flat=f.reals(key);
    if (flat.size()%static_cast<std::size_t>(n) ||
        flat.size()/static_cast<std::size_t>(n)<static_cast<std::size_t>(occupied))
        throw std::runtime_error("invalid MO coefficient dimensions");
    Matrix c(n,occupied);
    for (Eigen::Index j=0;j<occupied;++j)
        for (Eigen::Index i=0;i<n;++i)
            c(i,j)=flat[static_cast<std::size_t>(j*n+i)];
    return c;
}
template<class F> double seconds(int repeats, F&& fn) {
    fn();
    double best=std::numeric_limits<double>::infinity();
    for (int i=0;i<repeats;++i) {
        const auto t=Clock::now();
        fn();
        best=std::min(best,std::chrono::duration<double>(Clock::now()-t).count());
    }
    return best;
}
void check_close(const Eigen::VectorXd& p, const Eigen::VectorXd& mo,
                 double& max_abs, double& max_rel) {
    for (Eigen::Index i=0;i<p.size();++i) {
        if (!std::isfinite(p[i]) || !std::isfinite(mo[i]))
            throw std::runtime_error("nonfinite density");
        const double delta=std::abs(p[i]-mo[i]);
        max_abs=std::max(max_abs,delta);
        max_rel=std::max(max_rel,delta/std::max(1e-10,std::abs(p[i])));
    }
}
}
int main(int argc,char** argv) {
    if (argc<2 || argc>5) {
        std::cerr << "Usage: kumo-backend-comparison FILE.fchk[.gz] [points=4096] [repeats=5] [block=128]\n";
        return 2;
    }
    try {
        const std::size_t requested=argc>2?std::stoull(argv[2]):4096;
        const int repeats=argc>3?std::stoi(argv[3]):5;
        const std::size_t block=argc>4?std::stoull(argv[4]):128;
        if (requested<8 || requested>1000000 || repeats<1 || repeats>50 || block<1 || block>65536)
            throw std::invalid_argument("invalid benchmark parameters");
        const FormattedCheckpoint f(argv[1]);
        const chem::kumo::BasisEvaluator basis{chem::kumo::BasisSet{f}};
        const Eigen::Index n=static_cast<Eigen::Index>(basis.basis_size());
        const Eigen::Index na=f.integer("Number of alpha electrons");
        const Eigen::Index nb=f.integer("Number of beta electrons");
        if (na<0 || nb<0 || na>n || nb>n) throw std::runtime_error("invalid occupations");
        const auto p=unpack(f.reals("Total SCF Density"),n);
        const auto ca=coefficients(f,"Alpha MO coefficients",n,na);
        const bool has_beta=f.contains("Beta MO coefficients");
        if (!has_beta && na!=nb) throw std::runtime_error("missing beta MOs for open shell");
        const Matrix cb=has_beta?coefficients(f,"Beta MO coefficients",n,nb):ca;
        // Reconstruct the AO density matrix using the same integer occupations
        // assumed by the MO benchmark. This check is independent of the grid.
        Matrix reconstructed=ca*ca.transpose();
        if (has_beta) reconstructed.noalias() += cb*cb.transpose();
        else reconstructed *= 2.0;
        const double p_max_abs=(p-reconstructed).cwiseAbs().maxCoeff();
        const double p_fro_rel=(p-reconstructed).norm()/std::max(p.norm(),1e-300);
        const double p_max_rel=p_max_abs/std::max(p.cwiseAbs().maxCoeff(),1e-300);
        // Fixed, deterministic Cartesian points, spanning a box around the molecular geometry.
        const auto& xyz=f.reals("Current cartesian coordinates");
        if (xyz.size()<3 || xyz.size()%3) throw std::runtime_error("invalid coordinates");
        Eigen::Vector3d lo(xyz[0],xyz[1],xyz[2]), hi=lo;
        for (std::size_t i=0;i<xyz.size();i+=3) {
            const Eigen::Vector3d r(xyz[i],xyz[i+1],xyz[i+2]);
            lo=lo.cwiseMin(r); hi=hi.cwiseMax(r);
        }
        lo.array()-=2.5; hi.array()+=2.5;
        const std::size_t side=static_cast<std::size_t>(std::ceil(std::cbrt(static_cast<double>(requested))));
        const std::size_t points=side*side*side;
        Eigen::Matrix<double,3,Eigen::Dynamic> coords(3,static_cast<Eigen::Index>(points));
        for (std::size_t i=0;i<points;++i) {
            const std::size_t ix=i/(side*side), iy=(i/side)%side, iz=i%side;
            coords.col(static_cast<Eigen::Index>(i))=lo+(hi-lo).cwiseProduct(Eigen::Vector3d(
                double(ix)/double(side-1),double(iy)/double(side-1),double(iz)/double(side-1)));
        }
        // AO evaluation is measured independently. All kernel-only timings reuse precomputed AO values.
        Matrix phi(n,static_cast<Eigen::Index>(points));
        const double t_ao=seconds(repeats,[&] {
            for (std::size_t start=0;start<points;start+=block) {
                const auto b=static_cast<Eigen::Index>(std::min(block,points-start));
                Matrix tmp(n,b);
                basis.evaluate_batch(coords.middleCols(static_cast<Eigen::Index>(start),b),tmp);
                phi.middleCols(static_cast<Eigen::Index>(start),b)=tmp;
            }
        });
        Eigen::VectorXd rho_p(static_cast<Eigen::Index>(points));
        Eigen::VectorXd rho_mo(static_cast<Eigen::Index>(points));
        volatile double guard=0;
        const auto kernel_p=[&] {
            for (std::size_t start=0;start<points;start+=block) {
                const auto b=static_cast<Eigen::Index>(std::min(block,points-start));
                const auto x=phi.middleCols(static_cast<Eigen::Index>(start),b);
                const Matrix y=p*x;
                rho_p.segment(static_cast<Eigen::Index>(start),b)=(x.array()*y.array()).colwise().sum().transpose();
            }
            guard=guard+rho_p[static_cast<Eigen::Index>(points/2)];
        };
        const auto kernel_mo=[&] {
            for (std::size_t start=0;start<points;start+=block) {
                const auto b=static_cast<Eigen::Index>(std::min(block,points-start));
                const auto x=phi.middleCols(static_cast<Eigen::Index>(start),b);
                const Matrix a=ca.transpose()*x;
                if (has_beta) {
                    const Matrix beta=cb.transpose()*x;
                    rho_mo.segment(static_cast<Eigen::Index>(start),b)=(a.array().square().colwise().sum()+beta.array().square().colwise().sum()).transpose();
                } else {
                    rho_mo.segment(static_cast<Eigen::Index>(start),b)=2.0*a.array().square().colwise().sum().transpose().matrix();
                }
            }
            guard=guard+rho_mo[static_cast<Eigen::Index>(points/2)];
        };
        kernel_p(); kernel_mo();
        double max_abs=0,max_rel=0;
        check_close(rho_p,rho_mo,max_abs,max_rel);
        const double t_p=seconds(repeats,kernel_p);
        const double t_mo=seconds(repeats,kernel_mo);

        // True streaming paths: one reusable AO block and contraction buffers.
        // No full-grid AO cache is accessed inside these measurements.
        Matrix ao_block(n,static_cast<Eigen::Index>(block));
        Matrix p_block(n,static_cast<Eigen::Index>(block));
        Matrix alpha_block(na,static_cast<Eigen::Index>(block));
        Matrix beta_block(has_beta?nb:0,static_cast<Eigen::Index>(block));
        Eigen::VectorXd stream_p(static_cast<Eigen::Index>(points));
        Eigen::VectorXd stream_mo(static_cast<Eigen::Index>(points));
        const auto streaming_p=[&] {
            for (std::size_t start=0;start<points;start+=block) {
                const auto b=static_cast<Eigen::Index>(std::min(block,points-start));
                auto x=ao_block.leftCols(b);
                basis.evaluate_batch(coords.middleCols(static_cast<Eigen::Index>(start),b),x);
                p_block.leftCols(b).noalias()=p*x;
                stream_p.segment(static_cast<Eigen::Index>(start),b)=
                    (x.array()*p_block.leftCols(b).array()).colwise().sum().transpose();
            }
            guard=guard+stream_p[static_cast<Eigen::Index>(points/2)];
        };
        const auto streaming_mo=[&] {
            for (std::size_t start=0;start<points;start+=block) {
                const auto b=static_cast<Eigen::Index>(std::min(block,points-start));
                auto x=ao_block.leftCols(b);
                basis.evaluate_batch(coords.middleCols(static_cast<Eigen::Index>(start),b),x);
                alpha_block.leftCols(b).noalias()=ca.transpose()*x;
                if (has_beta) {
                    beta_block.leftCols(b).noalias()=cb.transpose()*x;
                    stream_mo.segment(static_cast<Eigen::Index>(start),b)=
                        (alpha_block.leftCols(b).array().square().colwise().sum()+
                         beta_block.leftCols(b).array().square().colwise().sum()).transpose();
                } else {
                    stream_mo.segment(static_cast<Eigen::Index>(start),b)=
                        2.0*alpha_block.leftCols(b).array().square().colwise().sum().transpose().matrix();
                }
            }
            guard=guard+stream_mo[static_cast<Eigen::Index>(points/2)];
        };
        streaming_p(); streaming_mo();
        double streaming_max_abs=0, streaming_max_rel=0;
        check_close(stream_p,stream_mo,streaming_max_abs,streaming_max_rel);
        const double t_stream_p=seconds(repeats,streaming_p);
        const double t_stream_mo=seconds(repeats,streaming_mo);
        const double flop_p=2.0*double(n)*double(n)*double(points);
        const double flop_mo=2.0*double(n)*double(na+nb)*double(points);
        std::cout<<std::setprecision(12)
                 <<"file="<<argv[1]<<" basis="<<n<<" alpha_occ="<<na<<" beta_occ="<<nb
                 <<" points="<<points<<" block="<<block<<" repeats="<<repeats<<"\n"
                 <<"ao_time_s="<<t_ao<<"\n"
                 <<"density_matrix_time_s="<<t_p<<" estimated_gflop="<<flop_p/1e9
                 <<" estimated_gflops="<<flop_p/t_p/1e9<<"\n"
                 <<"occupied_orbitals_time_s="<<t_mo<<" estimated_gflop="<<flop_mo/1e9
                 <<" estimated_gflops="<<flop_mo/t_mo/1e9<<"\n"
                 <<"kernel_speedup_mo_vs_p="<<t_p/t_mo<<" max_abs="<<max_abs
                 <<" max_rel_floor_1e-10="<<max_rel<<" checksum="<<rho_p.sum()<<"\n"
                 <<"density_matrix_reconstruction_max_abs="<<p_max_abs
                 <<" density_matrix_reconstruction_max_rel="<<p_max_rel
                 <<" density_matrix_reconstruction_fro_rel="<<p_fro_rel<<"\n"
                 <<"streaming_density_matrix_time_s="<<t_stream_p
                 <<" streaming_occupied_orbitals_time_s="<<t_stream_mo
                 <<" streaming_speedup_mo_vs_p="<<t_stream_p/t_stream_mo
                 <<" streaming_max_abs="<<streaming_max_abs
                 <<" streaming_max_rel_floor_1e-10="<<streaming_max_rel<<"\n";
        (void)guard;
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        return 1;
    }
}
