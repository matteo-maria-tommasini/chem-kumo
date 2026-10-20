#include "chem/kumo/local_density_overlap_symmetric.hpp"
#include "chem/kumo/local_density_overlap_symmetric_paired.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            for (Spin spin : {Spin::Total, Spin::Difference}) {
                const Eigen::MatrixXd p = molecule.density_matrix(spin);
                const auto t0 = std::chrono::steady_clock::now();
                const Eigen::MatrixXd reference = local_density_overlap_symmetric(basis, p);
                const auto t1 = std::chrono::steady_clock::now();
                const Eigen::MatrixXd paired = local_density_overlap_symmetric_paired(basis, p);
                const auto t2 = std::chrono::steady_clock::now();
                const double err = (reference - paired).cwiseAbs().maxCoeff();
                const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
                const double symmetry = (paired - paired.transpose()).cwiseAbs().maxCoeff();
                const double old_seconds = std::chrono::duration<double>(t1-t0).count();
                const double paired_seconds = std::chrono::duration<double>(t2-t1).count();
                std::cout << std::scientific << std::setprecision(12)
                          << "system=" << argv[arg]
                          << " spin=" << (spin == Spin::Total ? "total" : "difference")
                          << " nAO=" << p.rows()
                          << " max_abs=" << err << " max_scaled=" << err / scale
                          << " symmetry=" << symmetry
                          << " symmetric_seconds=" << old_seconds
                          << " paired_seconds=" << paired_seconds
                          << " speedup=" << old_seconds / paired_seconds << '\n';
                if (!paired.allFinite() || err > 1e-12 * scale || symmetry > 1e-12)
                    throw std::runtime_error("paired density overlap mismatch");
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
