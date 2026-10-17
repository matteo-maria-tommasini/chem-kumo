#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Clock = std::chrono::steady_clock;
using chem::kumo::FormattedCheckpoint;

Matrix unpack(const std::vector<double>& data, Eigen::Index n) {
    if (data.size() != static_cast<std::size_t>(n*(n+1)/2))
        throw std::runtime_error("invalid packed density matrix");
    Matrix p(n,n);
    std::size_t k=0;
    for (Eigen::Index i=0;i<n;++i)
        for (Eigen::Index j=0;j<=i;++j) { p(i,j)=data[k]; p(j,i)=data[k]; ++k; }
    return p;
}
Matrix occupied(const FormattedCheckpoint& f, const std::string& name,
                Eigen::Index n, Eigen::Index count) {
    const auto& raw=f.reals(name);
    if (raw.size()%static_cast<std::size_t>(n) ||
        raw.size()/static_cast<std::size_t>(n)<static_cast<std::size_t>(count))
        throw std::runtime_error("invalid MO coefficient dimensions");
    Matrix c(n,count);
    for (Eigen::Index j=0;j<count;++j)
        for (Eigen::Index i=0;i<n;++i)
            c(i,j)=raw[static_cast<std::size_t>(j*n+i)];
    return c;
}
template<class F> double time_best(int repeats, F&& f) {
    f();
    double best=std::numeric_limits<double>::infinity();
    for (int i=0;i<repeats;++i) {
        const auto start=Clock::now();
        f();
        best=std::min(best,std::chrono::duration<double>(Clock::now()-start).count());
    }
    return best;
}
// Components: rho, dx,dy,dz, tau, xx,yy,zz,xy,xz,yz.
constexpr int RHO=0, TAU=4, XX=5;
using Output = Eigen::Matrix<double,11,Eigen::Dynamic>;
struct AOBlock {
    std::array<Matrix,10> a; // value, dx,dy,dz, xx,yy,zz,xy,xz,yz
    AOBlock(Eigen::Index n, Eigen::Index b) {
        for (auto& m:a) m.resize(n,b);
    }
};
void evaluate_ao(const chem::kumo::BasisEvaluator& basis,
                 const Eigen::Matrix<double,3,Eigen::Dynamic>& coords,
                 std::size_t start, Eigen::Index b, int order, AOBlock& block) {
    for (Eigen::Index j=0;j<b;++j) {
        const auto d=basis.evaluate_derivatives(coords.col(static_cast<Eigen::Index>(start)+j));
        block.a[0].col(j)=d.values;
        for (int k=0;k<3;++k) block.a[1+k].col(j)=d.gradient.col(k);
        if (order==2)
            for (int k=0;k<6;++k) block.a[4+k].col(j)=d.hessian.col(k);
    }
}
void accumulate_p(const Matrix& p, const AOBlock& a, Eigen::Index b,
                  int mode, Output& out) {
    const auto x=a.a[0].leftCols(b);
    const Matrix y=p*x;
    out.row(RHO).head(b)=(x.array()*y.array()).colwise().sum();
    for (int k=0;k<3;++k)
        out.row(1+k).head(b)=2.0*(a.a[1+k].leftCols(b).array()*y.array()).colwise().sum();
    if (mode==0) return;
    std::array<Matrix,3> gy;
    for (int k=0;k<3;++k) gy[k]=p*a.a[1+k].leftCols(b);
    out.row(TAU).head(b).setZero();
    for (int k=0;k<3;++k)
        out.row(TAU).head(b).array() +=
            0.5*(a.a[1+k].leftCols(b).array()*gy[k].array()).colwise().sum().array();
    if (mode==1) return;
    constexpr int pairs[6][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
    for (int h=0;h<6;++h) {
        const auto [k,l]=std::pair{pairs[h][0],pairs[h][1]};
        out.row(XX+h).head(b)=2.0*(
            (a.a[4+h].leftCols(b).array()*y.array()).colwise().sum()+
            (a.a[1+k].leftCols(b).array()*gy[l].array()).colwise().sum());
    }
}
void accumulate_mo_spin(const Matrix& c, const AOBlock& a, Eigen::Index b,
                        int mode, Output& out, double weight) {
    const Matrix psi=c.transpose()*a.a[0].leftCols(b);
    out.row(RHO).head(b).array() += weight*psi.array().square().colwise().sum().array();
    std::array<Matrix,3> g;
    for (int k=0;k<3;++k) {
        g[k]=c.transpose()*a.a[1+k].leftCols(b);
        out.row(1+k).head(b).array() +=
            2.0*weight*(psi.array()*g[k].array()).colwise().sum().array();
    }
    if (mode==0) return;
    for (int k=0;k<3;++k)
        out.row(TAU).head(b).array() +=
            0.5*weight*g[k].array().square().colwise().sum().array();
    if (mode==1) return;
    constexpr int pairs[6][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
    for (int h=0;h<6;++h) {
        const auto d2=c.transpose()*a.a[4+h].leftCols(b);
        const auto [k,l]=std::pair{pairs[h][0],pairs[h][1]};
        out.row(XX+h).head(b).array() += 2.0*weight*(
            (psi.array()*d2.array()).colwise().sum().array()+
            (g[k].array()*g[l].array()).colwise().sum().array());
    }
}
}
int main(int argc, char** argv) {
    if (argc!=7 || std::string(argv[5])!="--mode") {
        std::cerr << "Usage: kumo-derivative-backend-comparison FILE points repeats block --mode density-gradient|elf-ingredients|density-hessian\n";
        return 2;
    }
    try {
        const std::size_t requested=std::stoull(argv[2]);
        const int repeats=std::stoi(argv[3]);
        const std::size_t block_size=std::stoull(argv[4]);
        const std::string mode_name=argv[6];
        const int mode=mode_name=="density-gradient"?0:
                       mode_name=="elf-ingredients"?1:
                       mode_name=="density-hessian"?2:-1;
        if (mode<0 || requested<8 || requested>100000 || repeats<1 || repeats>50 ||
            block_size<1 || block_size>4096)
            throw std::runtime_error("invalid benchmark parameters");
        const FormattedCheckpoint f(argv[1]);
        const chem::kumo::BasisEvaluator basis{chem::kumo::BasisSet{f}};
        const Eigen::Index n=static_cast<Eigen::Index>(basis.basis_size());
        const Eigen::Index na=f.integer("Number of alpha electrons");
        const Eigen::Index nb=f.integer("Number of beta electrons");
        if (na<0 || nb<0 || na>n || nb>n) throw std::runtime_error("invalid electron counts");
        const Matrix p=unpack(f.reals("Total SCF Density"),n);
        const Matrix ca=occupied(f,"Alpha MO coefficients",n,na);
        const bool unrestricted=f.contains("Beta MO coefficients");
        if (!unrestricted && na!=nb) throw std::runtime_error("open shell without beta coefficients");
        const Matrix cb=unrestricted?occupied(f,"Beta MO coefficients",n,nb):Matrix();
        const auto& xyz=f.reals("Current cartesian coordinates");
        if (xyz.empty() || xyz.size()%3) throw std::runtime_error("invalid geometry");
        Eigen::Vector3d lo(xyz[0],xyz[1],xyz[2]), hi=lo;
        for (std::size_t i=0;i<xyz.size();i+=3) {
            Eigen::Vector3d r(xyz[i],xyz[i+1],xyz[i+2]);
            lo=lo.cwiseMin(r); hi=hi.cwiseMax(r);
        }
        lo.array()-=2.5; hi.array()+=2.5;
        const std::size_t side=static_cast<std::size_t>(std::ceil(std::cbrt(double(requested))));
        const std::size_t points=side*side*side;
        Eigen::Matrix<double,3,Eigen::Dynamic> coords(3,static_cast<Eigen::Index>(points));
        for (std::size_t i=0;i<points;++i) {
            const std::size_t ix=i/(side*side),iy=(i/side)%side,iz=i%side;
            coords.col(static_cast<Eigen::Index>(i))=lo+(hi-lo).cwiseProduct(Eigen::Vector3d(
                double(ix)/double(side-1),double(iy)/double(side-1),double(iz)/double(side-1)));
        }
        AOBlock ao(n,static_cast<Eigen::Index>(block_size));
        Output result_p(11,static_cast<Eigen::Index>(points));
        Output result_mo(11,static_cast<Eigen::Index>(points));
        Output local(11,static_cast<Eigen::Index>(block_size));
        const int order=mode==2?2:1;
        auto run=[&](bool use_mo, bool evaluate, double& sink) {
            auto& output=use_mo?result_mo:result_p;
            for (std::size_t start=0;start<points;start+=block_size) {
                const auto b=static_cast<Eigen::Index>(std::min(block_size,points-start));
                if (evaluate) evaluate_ao(basis,coords,start,b,order,ao);
                local.leftCols(b).setZero();
                if (use_mo) {
                    accumulate_mo_spin(ca,ao,b,mode,local,unrestricted?1.0:2.0);
                    if (unrestricted) accumulate_mo_spin(cb,ao,b,mode,local,1.0);
                } else accumulate_p(p,ao,b,mode,local);
                output.middleCols(static_cast<Eigen::Index>(start),b)=local.leftCols(b);
            }
            sink += output(RHO,static_cast<Eigen::Index>(points/2));
        };
        double sink=0.0;
        run(false,true,sink); run(true,true,sink);
        const int rows=mode==0?4:mode==1?5:11;
        double max_abs=0.0,max_scaled=0.0;
        for (int k=0;k<rows;++k)
            for (Eigen::Index j=0;j<static_cast<Eigen::Index>(points);++j) {
                const double a=result_p(k,j),b=result_mo(k,j);
                if (!std::isfinite(a) || !std::isfinite(b)) throw std::runtime_error("nonfinite output");
                max_abs=std::max(max_abs,std::abs(a-b));
                max_scaled=std::max(max_scaled,std::abs(a-b)/std::max(1e-10,std::abs(a)));
            }
        // Stage 01D: property-resolved diagnostics on the same completed grids.
        // Components use the ordering declared by Output above.
        constexpr std::array<const char*,11> names={
            "rho", "grad_x", "grad_y", "grad_z", "tau",
            "hess_xx", "hess_yy", "hess_zz", "hess_xy", "hess_xz", "hess_yz"};
        std::array<double,11> component_abs{}, component_scaled{};
        std::array<std::size_t,11> component_worst{};
        for (int k=0;k<rows;++k) {
            for (Eigen::Index j=0;j<static_cast<Eigen::Index>(points);++j) {
                const double a=result_p(k,j), b=result_mo(k,j);
                const double err=std::abs(a-b);
                if (err>component_abs[k]) {
                    component_abs[k]=err;
                    component_worst[k]=static_cast<std::size_t>(j);
                }
                component_scaled[k]=std::max(component_scaled[k],
                    err/std::max(1e-10,std::abs(a)));
            }
        }
        // An independent eigenvalue calculation for each Hessian representation.
        // Weyl's inequality bounds each ordered eigenvalue error by the
        // spectral norm of the matrix difference.
        double eigen_max_abs=0.0, eigen_max_bound_excess=0.0;
        double laplacian_max_abs=0.0, eigen_trace_max_abs=0.0;
        std::size_t lambda2_sign_ambiguous=0, lambda2_sign_stable_mismatch=0;
        std::array<std::size_t,3> low_density_counts{};
        constexpr std::array<double,3> low_density_cutoffs={1e-4,1e-8,1e-12};
        for (Eigen::Index j=0;j<static_cast<Eigen::Index>(points);++j) {
            const double rho=result_p(RHO,j);
            for (std::size_t t=0;t<low_density_cutoffs.size();++t)
                if (std::abs(rho)<low_density_cutoffs[t]) ++low_density_counts[t];
            if (mode!=2) continue;
            const auto make_hessian=[&](const Output& o) {
                Eigen::Matrix3d h;
                h << o(5,j),o(8,j),o(9,j),
                     o(8,j),o(6,j),o(10,j),
                     o(9,j),o(10,j),o(7,j);
                return h;
            };
            const Eigen::Matrix3d hp=make_hessian(result_p);
            const Eigen::Matrix3d hm=make_hessian(result_mo);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> ep(hp), em(hm);
            if (ep.info()!=Eigen::Success || em.info()!=Eigen::Success)
                throw std::runtime_error("Hessian diagonalization failed");
            const double bound=(hp-hm).norm(); // Frobenius >= spectral norm
            for (int k=0;k<3;++k) {
                const double err=std::abs(ep.eigenvalues()[k]-em.eigenvalues()[k]);
                eigen_max_abs=std::max(eigen_max_abs,err);
                eigen_max_bound_excess=std::max(eigen_max_bound_excess,err-bound);
            }
            const double lp=hp.trace(), lm=hm.trace();
            laplacian_max_abs=std::max(laplacian_max_abs,std::abs(lp-lm));
            eigen_trace_max_abs=std::max(eigen_trace_max_abs,
                std::max(std::abs(ep.eigenvalues().sum()-lp),
                         std::abs(em.eigenvalues().sum()-lm)));
            const double l2p=ep.eigenvalues()[1], l2m=em.eigenvalues()[1];
            if (std::abs(l2p)<=bound || std::abs(l2m)<=bound)
                ++lambda2_sign_ambiguous;
            else if (std::signbit(l2p)!=std::signbit(l2m))
                ++lambda2_sign_stable_mismatch;
        }
        if (eigen_max_bound_excess>1e-10 || lambda2_sign_stable_mismatch!=0)
            throw std::runtime_error("Hessian eigenvalue stability invariant failed");
        const double t_ao=time_best(repeats,[&] {
            for (std::size_t start=0;start<points;start+=block_size) {
                const auto b=static_cast<Eigen::Index>(std::min(block_size,points-start));
                evaluate_ao(basis,coords,start,b,order,ao);
            }
            sink+=ao.a[0](0,0);
        });
        const double t_p=time_best(repeats,[&] { run(false,true,sink); });
        const double t_mo=time_best(repeats,[&] { run(true,true,sink); });
        // Timings include pointwise AO derivatives and fused streaming contraction.
        std::cout<<std::setprecision(12)<<"file="<<argv[1]<<" mode="<<mode_name
                 <<" basis="<<n<<" alpha_occ="<<na<<" beta_occ="<<nb
                 <<" points="<<points<<" block="<<block_size<<" repeats="<<repeats<<'\n'
                 <<"ao_derivatives_time_s="<<t_ao<<'\n'
                 <<"streaming_density_matrix_time_s="<<t_p<<'\n'
                 <<"streaming_occupied_orbitals_time_s="<<t_mo<<'\n'
                 <<"streaming_speedup_mo_vs_p="<<t_p/t_mo<<'\n'
                 <<"max_abs="<<max_abs<<" max_rel_floor_1e-10="<<max_scaled
                 <<" checksum="<<result_p.row(RHO).sum()<<'\n';
        for (int k=0;k<rows;++k)
            std::cout<<"component="<<names[k]
                     <<" max_abs="<<component_abs[k]
                     <<" max_rel_floor_1e-10="<<component_scaled[k]
                     <<" worst_point="<<component_worst[k]<<'\n';
        for (std::size_t t=0;t<low_density_cutoffs.size();++t)
            std::cout<<"low_density_threshold="<<low_density_cutoffs[t]
                     <<" points="<<low_density_counts[t]<<'\n';
        if (mode==2)
            std::cout<<"hessian_eigen_max_abs="<<eigen_max_abs
                     <<" hessian_eigen_bound_excess="<<eigen_max_bound_excess
                     <<" laplacian_max_abs="<<laplacian_max_abs
                     <<" eigen_trace_max_abs="<<eigen_trace_max_abs
                     <<" lambda2_sign_ambiguous="<<lambda2_sign_ambiguous
                     <<" lambda2_sign_stable_mismatch="<<lambda2_sign_stable_mismatch<<'\n';
        (void)sink;
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        return 1;
    }
}
