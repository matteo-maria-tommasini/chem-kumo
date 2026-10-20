#include "chem/kumo/local_density_overlap.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

static double maxabs(const Eigen::MatrixXd& x) {
    return x.size() ? x.cwiseAbs().maxCoeff() : 0.0;
}

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto pa = molecule.density_matrix(Spin::Alpha);
            const auto pb = molecule.density_matrix(Spin::Beta);
            const auto pt = molecule.density_matrix(Spin::Total);
            const auto ps = molecule.density_matrix(Spin::Difference);
            const auto wa = local_density_overlap(basis, pa);
            const auto wb = local_density_overlap(basis, pb);
            const auto wt = local_density_overlap(basis, pt);
            const auto ws = local_density_overlap(basis, ps);
            const double linear = maxabs(wt - wa - wb);
            const double spin_linear = maxabs(ws - wa + wb);
            const double symmetry = maxabs(wt - wt.transpose());
            require(linear < 1e-9 && spin_linear < 1e-9, "spin linearity failure");
            require(symmetry < 1e-12 && wt.allFinite(), "symmetry/finite failure");
            // Positivity of the weighted overlap for a nonnegative SCF density.
            for (Eigen::Index k = 0; k < std::min<Eigen::Index>(wt.rows(), 12); ++k) {
                Eigen::VectorXd v(wt.rows());
                for (Eigen::Index i = 0; i < v.size(); ++i)
                    v[i] = std::sin(0.37 * (i + 1) * (k + 2));
                require(v.dot(wt * v) > -1e-8 * v.squaredNorm(), "PSD failure");
            }
            const double rho2 = (pt.cwiseProduct(wt)).sum();
            require(std::isfinite(rho2) && rho2 > 0.0, "invalid integral rho squared");
            // A nonsymmetric input must be rejected rather than silently altered.
            Eigen::MatrixXd bad = pt;
            if (bad.rows() > 1) {
                bad(0,1) += 1e-2;
                bool rejected = false;
                try { (void)local_density_overlap(basis, bad); }
                catch (const std::invalid_argument&) { rejected = true; }
                require(rejected, "asymmetric input accepted");
            }
            std::cout << argv[arg] << " nAO=" << wt.rows()
                      << " linearity=" << std::scientific << std::setprecision(12) << linear
                      << " spin_linearity=" << spin_linear
                      << " symmetry=" << symmetry
                      << " integral_rho2=" << rho2 << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
