#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void check_close(double a, double b, const char* label) {
    if (!std::isfinite(a) || !std::isfinite(b) ||
        std::abs(a-b) > 2.e-5 + 2.e-5*std::max(std::abs(a),std::abs(b)))
        throw std::runtime_error(label);
}
void compare(const chem::kumo::DensityDerivatives& a,
             const chem::kumo::DensityDerivatives& b) {
    check_close(a.value,b.value,"MO/DM density mismatch");
    check_close(a.tau,b.tau,"MO/DM tau mismatch");
    check_close(a.laplacian,b.laplacian,"MO/DM Laplacian mismatch");
    for (int i=0;i<3;++i) {
        check_close(a.gradient[i],b.gradient[i],"MO/DM gradient mismatch");
        for (int j=0;j<3;++j)
            check_close(a.hessian(i,j),b.hessian(i,j),"MO/DM Hessian mismatch");
    }
}
chem::kumo::DensityDerivatives from_mos(const chem::kumo::Molecule& mol,
                                         const Eigen::Vector3d& point,
                                         chem::kumo::Spin spin) {
    chem::kumo::DensityDerivatives out;
    const int occupied = mol.occupied_count(spin);
    for (int i=1;i<=occupied;++i) {
        const auto mo=mol.orbital_derivatives(point,static_cast<std::size_t>(i),spin);
        out.value += mo.value*mo.value;
        out.gradient += 2.0*mo.value*mo.gradient;
        out.hessian += 2.0*(mo.gradient*mo.gradient.transpose()+mo.value*mo.hessian);
        out.tau += 0.5*mo.gradient.squaredNorm();
    }
    out.laplacian=out.hessian.trace();
    return out;
}
}
int main(int argc, char** argv) {
    if (argc!=3) return 2;
    for (int f=1;f<argc;++f) {
        const auto mol=chem::kumo::Molecule::from_fchk(argv[f]);
        Eigen::Matrix<double,3,Eigen::Dynamic> points(3,13);
        for (int j=0;j<points.cols();++j)
            points.col(j)=Eigen::Vector3d(0.31+0.21*j,-0.47+0.17*j,0.39-0.11*j);
        for (auto spin : {chem::kumo::Spin::Alpha,chem::kumo::Spin::Beta}) {
            const auto block=mol.density_derivatives_blocked(points,spin,5);
            const auto block1=mol.density_derivatives_blocked(points,spin,1);
            if (block.size()!=static_cast<std::size_t>(points.cols()))
                throw std::runtime_error("wrong block output size");
            for (int j=0;j<points.cols();++j) {
                const auto r=points.col(j).eval();
                const auto reference=mol.density_derivatives(r,spin);
                compare(block[static_cast<std::size_t>(j)],reference);
                compare(block1[static_cast<std::size_t>(j)],reference);
                compare(from_mos(mol,r,spin),reference);
            }
        }
        for (auto spin : {chem::kumo::Spin::Total,chem::kumo::Spin::Difference}) {
            const auto block=mol.density_derivatives_blocked(points,spin,7);
            for (int j=0;j<points.cols();++j)
                compare(block[static_cast<std::size_t>(j)],
                        mol.density_derivatives(points.col(j),spin));
        }
        bool rejected=false;
        try { (void)mol.density_derivatives_blocked(points,chem::kumo::Spin::Total,0); }
        catch (const std::invalid_argument&) { rejected=true; }
        if (!rejected) throw std::runtime_error("zero block size accepted");
        std::cout << argv[f] << " MO/DM and blocked derivatives OK\n";
    }
}
