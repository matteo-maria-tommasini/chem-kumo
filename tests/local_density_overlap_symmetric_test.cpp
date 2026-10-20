#include "chem/kumo/local_density_overlap.hpp"
#include "chem/kumo/local_density_overlap_symmetric.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

static double maxabs(const Eigen::MatrixXd& a) {
    return a.size() ? a.cwiseAbs().maxCoeff() : 0.0;
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            for (Spin spin : {Spin::Total, Spin::Difference}) {
                const Eigen::MatrixXd p = molecule.density_matrix(spin);
                const auto t0 = std::chrono::steady_clock::now();
                const Eigen::MatrixXd reference = local_density_overlap(basis, p);
                const auto t1 = std::chrono::steady_clock::now();
                const Eigen::MatrixXd optimized = local_density_overlap_symmetric(basis, p);
                const auto t2 = std::chrono::steady_clock::now();
                const double err = maxabs(reference - optimized);
                const double scale = std::max(1.0, maxabs(reference));
                const double symmetry = maxabs(optimized - optimized.transpose());
                const double direct_seconds = std::chrono::duration<double>(t1-t0).count();
                const double symmetric_seconds = std::chrono::duration<double>(t2-t1).count();
                std::cout << std::scientific << std::setprecision(12)
                          << argv[arg] << " spin=" << (spin == Spin::Total ? "total" : "difference")
                          << " nAO=" << p.rows() << " max_abs=" << err
                          << " max_scaled=" << err/scale << " symmetry=" << symmetry
                          << " direct_seconds=" << direct_seconds
                          << " symmetric_seconds=" << symmetric_seconds
                          << " speedup=" << direct_seconds/symmetric_seconds << '\n';
                if (!optimized.allFinite() || err > 1e-11 * scale || symmetry > 1e-12)
                    throw std::runtime_error("symmetry-reduced density overlap mismatch");
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
