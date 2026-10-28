#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
constexpr double au_to_debye = 2.541746473;

int validate(const std::string& path) {
    using namespace chem::kumo;
    const FormattedCheckpoint fchk(path);
    if (!fchk.contains("Dipole Moment"))
        throw std::runtime_error("Gaussian Dipole Moment record absent: " + path);
    const auto& reference = fchk.reals("Dipole Moment");
    if (reference.size() != 3)
        throw std::runtime_error("Gaussian Dipole Moment must have 3 Cartesian components");

    const Molecule mol = Molecule::from_fchk(path);
    const BasisEvaluator ao{BasisSet(fchk)};
    const CartesianAOMoments moments = ao.cartesian_moments();
    const Eigen::MatrixXd P = mol.density_matrix(Spin::Total);
    if (P.rows() != moments.overlap.rows() || P.cols() != moments.overlap.cols())
        throw std::runtime_error("P/AO basis dimension mismatch");

    Eigen::Vector3d nuclear = Eigen::Vector3d::Zero();
    double nuclear_charge = 0.0;
    for (const auto& atom : mol.atoms()) {
        nuclear += atom.nuclear_charge * atom.position;
        nuclear_charge += atom.nuclear_charge;
    }
    Eigen::Vector3d electron = Eigen::Vector3d::Zero();
    for (int k = 0; k < 3; ++k)
        electron[k] = (P.cwiseProduct(moments.first[static_cast<std::size_t>(k)])).sum();
    const Eigen::Vector3d dipole = nuclear - electron;
    const Eigen::Vector3d gaussian(reference[0], reference[1], reference[2]);
    const double electrons = (P.cwiseProduct(moments.overlap)).sum();
    const double charge = nuclear_charge - electrons;
    const double abs_error = (dipole - gaussian).cwiseAbs().maxCoeff();
    const double tolerance = 2e-5; // a.u.; Gaussian FCHK precision and AO contractions

    std::cout << std::setprecision(12)
              << "file=" << path << " N=" << electrons << " Q=" << charge << '\n'
              << "  nuclear_au=" << nuclear.transpose() << '\n'
              << "  electronic_au=" << electron.transpose() << '\n'
              << "  kumo_au=" << dipole.transpose() << '\n'
              << "  gaussian_au=" << gaussian.transpose() << '\n'
              << "  delta_au=" << (dipole - gaussian).transpose() << '\n'
              << "  max_abs_error_au=" << abs_error
              << "  max_abs_error_D=" << abs_error * au_to_debye << '\n';
    if (!dipole.allFinite() || abs_error > tolerance)
        throw std::runtime_error("dipole differs from Gaussian: " + path);
    return 0;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-gaussian-dipole-reference-test file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        for (int i = 1; i < argc; ++i) validate(argv[i]);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
