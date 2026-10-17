#include "chem/kumo/density.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace chem::kumo::detail {
struct AoMoTestAccess {
    static auto run(const Molecule& m,
        const Eigen::Ref<const Eigen::Matrix<double,3,Eigen::Dynamic>>& points,
        const std::vector<Spin>& components, std::size_t block, int order) {
        return m.density_derivatives_mo(points, components, block, order);
    }
    static auto matrices(const Molecule& m) {
        const auto& cb = m.unrestricted_ ? m.beta_mos_ : m.alpha_mos_;
        const Eigen::MatrixXd pa = m.alpha_mos_.leftCols(m.nalpha_) *
                                   m.alpha_mos_.leftCols(m.nalpha_).transpose();
        const Eigen::MatrixXd pb = cb.leftCols(m.nbeta_) * cb.leftCols(m.nbeta_).transpose();
        std::map<Spin,Eigen::MatrixXd> reconstructed;
        reconstructed.emplace(Spin::Alpha,pa);
        reconstructed.emplace(Spin::Beta,pb);
        reconstructed.emplace(Spin::Total,pa+pb);
        reconstructed.emplace(Spin::Difference,pa-pb);
        std::map<Spin,Eigen::MatrixXd> original;
        original.emplace(Spin::Total,m.total_);
        original.emplace(Spin::Difference,m.spin_);
        original.emplace(Spin::Alpha,0.5*(m.total_+m.spin_));
        original.emplace(Spin::Beta,0.5*(m.total_-m.spin_));
        return std::make_pair(std::move(reconstructed),std::move(original));
    }
    static auto reference(const Molecule& m, const Eigen::Vector3d& point,
                          const Eigen::MatrixXd& p) {
        const auto ao = m.evaluator_.evaluate_derivatives(point,2);
        const auto y = p * ao.values;
        const auto pg = p * ao.gradient;
        constexpr int index[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
        DensityDerivatives d;
        d.value = ao.values.dot(y);
        for (int a=0;a<3;++a) {
            d.gradient[a] = 2.0*ao.gradient.col(a).dot(y);
            d.tau += 0.5*ao.gradient.col(a).dot(pg.col(a));
            for (int b=a;b<3;++b) {
                d.hessian(a,b)=d.hessian(b,a)=2.0*(
                    ao.hessian.col(index[a][b]).dot(y) +
                    ao.gradient.col(a).dot(pg.col(b)));
            }
        }
        d.laplacian = d.hessian.trace();
        return d;
    }
};
}

using namespace chem::kumo;
namespace {
struct ErrorRecord {
    double max_abs=0.0;
    double max_scaled=0.0;
    double left=0.0;
    double right=0.0;
    std::size_t point=0;
    std::string component;
};
void accumulate(ErrorRecord& rec, double left, double right,
                std::size_t point, const std::string& component) {
    if (!std::isfinite(left) || !std::isfinite(right))
        throw std::runtime_error("nonfinite diagnostic value");
    const double diff=std::abs(left-right);
    if (diff>=rec.max_abs) {
        rec.max_abs=diff; rec.left=left; rec.right=right;
        rec.point=point; rec.component=component;
    }
    rec.max_scaled=std::max(rec.max_scaled,
        diff/std::max({1.0,std::abs(left),std::abs(right)}));
}
struct Comparison {
    ErrorRecord rho,grad,hessian,laplacian,tau;
    void add(const DensityDerivatives& a, const DensityDerivatives& b,
             std::size_t j) {
        accumulate(rho,a.value,b.value,j,"rho");
        for (int k=0;k<3;++k)
            accumulate(grad,a.gradient[k],b.gradient[k],j,
                       std::string("grad_")+"xyz"[k]);
        constexpr int pairs[6][2]={{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
        constexpr const char* names[6]={"xx","yy","zz","xy","xz","yz"};
        for (int k=0;k<6;++k)
            accumulate(hessian,a.hessian(pairs[k][0],pairs[k][1]),
                       b.hessian(pairs[k][0],pairs[k][1]),j,names[k]);
        accumulate(laplacian,a.laplacian,b.laplacian,j,"laplacian");
        accumulate(tau,a.tau,b.tau,j,"tau");
    }
    void print(const char* label) const {
        std::cout << "  " << label << '\n';
        for (const auto& [name,rec]: std::array<std::pair<const char*,const ErrorRecord*>,5>{{
            {"rho",&rho},{"gradient",&grad},{"hessian",&hessian},
            {"laplacian",&laplacian},{"tau",&tau}}}) {
            std::cout << "    " << name << " max_abs=" << rec->max_abs
                      << " max_scaled=" << rec->max_scaled
                      << " point=" << rec->point << " component=" << rec->component
                      << " lhs=" << rec->left << " rhs=" << rec->right << '\n';
        }
    }
};
const char* spin_name(Spin s) {
    switch(s) {
        case Spin::Total:return "total";
        case Spin::Alpha:return "alpha";
        case Spin::Beta:return "beta";
        case Spin::Difference:return "spin";
    }
    return "unknown";
}
}
int main(int argc,char** argv) {
    if (argc<2) {
        std::cerr << "usage: kumo-ao-mo-diagnostic file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(16);
        const std::vector<Spin> spins{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference};
        Eigen::Matrix<double,3,Eigen::Dynamic> points(3,19);
        for (int j=0;j<19;++j)
            points.col(j)=Eigen::Vector3d(0.17*j-1.2,0.09*j-0.8,0.13*j+0.2);
        for (int arg=1;arg<argc;++arg) {
            const auto m=Molecule::from_fchk(argv[arg]);
            if (!m.electronic_structure_info().ao_mo_eligible())
                throw std::runtime_error(std::string("AO-MO ineligible: ")+argv[arg]);
            const auto [reconstructed,original]=detail::AoMoTestAccess::matrices(m);
            std::cout << "file=" << argv[arg] << " basis=" << m.basis_size() << '\n';
            for (Spin spin:spins) {
                const auto& pm=reconstructed.at(spin);
                const auto& pf=original.at(spin);
                std::cout << "spin=" << spin_name(spin)
                          << " matrix_frobenius=" << (pm-pf).norm() << '\n';
            }
            for (std::size_t block:{std::size_t{1},std::size_t{7},std::size_t{32}}) {
                const auto mo=detail::AoMoTestAccess::run(m,points,spins,block,2);
                std::cout << "block=" << block << '\n';
                for (Spin spin:spins) {
                    Comparison algorithm,representation,total;
                    for (std::size_t j=0;j<static_cast<std::size_t>(points.cols());++j) {
                        const auto rec=detail::AoMoTestAccess::reference(
                            m,points.col(static_cast<Eigen::Index>(j)),reconstructed.at(spin));
                        const auto fchk=detail::AoMoTestAccess::reference(
                            m,points.col(static_cast<Eigen::Index>(j)),original.at(spin));
                        const auto& orbital=mo.at(spin).at(j);
                        algorithm.add(orbital,rec,j);
                        representation.add(rec,fchk,j);
                        total.add(orbital,fchk,j);
                    }
                    std::cout << " spin=" << spin_name(spin) << '\n';
                    algorithm.print("algorithm: MO vs AO-P(P_MO)");
                    representation.print("representation: AO-P(P_MO) vs AO-P(P_FCHK)");
                    total.print("observed: MO vs AO-P(P_FCHK)");
                }
            }
        }
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "diagnostic error: " << e.what() << '\n';
        return 1;
    }
}
