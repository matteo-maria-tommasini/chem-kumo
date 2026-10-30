#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Core>
#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
int check(const std::string& path) {
    using namespace chem::kumo;
    const FormattedCheckpoint fchk(path);
    if (!fchk.contains("Quadrupole Moment"))
        throw std::runtime_error("Gaussian Quadrupole Moment absent: " + path);
    const auto& r = fchk.reals("Quadrupole Moment");
    if (r.size() != 6) throw std::runtime_error("Quadrupole Moment must have six components");
    const Molecule mol = Molecule::from_fchk(path);
    const BasisEvaluator ao{BasisSet(fchk)};
    const auto moments = ao.cartesian_moments();
    const Eigen::MatrixXd P = mol.density_matrix(Spin::Total);
    if (P.rows() != moments.overlap.rows() || P.cols() != moments.overlap.cols())
        throw std::runtime_error("P / AO dimension mismatch");
    constexpr int iidx[6] = {0,1,2,0,0,1};
    constexpr int jidx[6] = {0,1,2,1,2,2};
    std::array<double,6> primitive{};
    for (int k=0; k<6; ++k) {
        double nuclear = 0;
        for (const auto& atom : mol.atoms())
            nuclear += atom.nuclear_charge * atom.position[iidx[k]] * atom.position[jidx[k]];
        primitive[k] = nuclear - (P.cwiseProduct(moments.second[k])).sum();
    }
    const double tr = primitive[0]+primitive[1]+primitive[2];
    std::array<double,6> gaussian{};
    std::array<double,6> conventional{};
    std::array<double,6> buckingham{};
    double maxerr=0;
    double maxerr_buck=0;
    for (int k=0;k<6;++k) {
        gaussian[k]=r[k];
        conventional[k]=primitive[k] - (k<3 ? tr/3.0 : 0.0);
        buckingham[k]=1.5*conventional[k];
        maxerr=std::max(maxerr,std::abs(conventional[k]-gaussian[k]));
        maxerr_buck=std::max(maxerr_buck,std::abs(buckingham[k]-gaussian[k]));
    }
    std::cout << std::setprecision(12) << "file=" << path << '\n';
    for(int k=0;k<6;++k)
        std::cout << "  " << std::array<const char*,6>{"xx","yy","zz","xy","xz","yz"}[k]
                  << " primitive=" << primitive[k] << " kumo_traceless=" << conventional[k]
                  << " gaussian=" << gaussian[k] << " delta=" << conventional[k]-gaussian[k] << '\n';
    std::cout << "  reference_trace=" << (gaussian[0]+gaussian[1]+gaussian[2])
              << " max_abs_error_au=" << maxerr
              << " buckingham_max_abs_error_au=" << maxerr_buck << '\n';
    if (!std::isfinite(maxerr) || maxerr > 2e-5)
        throw std::runtime_error("quadrupole differs from Gaussian: " + path);
    return 0;
}
}
int main(int argc,char** argv) {
    if(argc<2) {std::cerr << "usage: kumo-gaussian-quadrupole-reference-test file.fchk[.gz] ...\n"; return 2;}
    try { for(int i=1;i<argc;++i) check(argv[i]); }
    catch(const std::exception& e) {std::cerr << "ERROR: " << e.what() << '\n'; return 1;}
    return 0;
}
