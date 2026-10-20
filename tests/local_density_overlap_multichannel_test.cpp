#include "chem/kumo/local_density_overlap_multichannel.hpp"
#include "chem/kumo/local_density_overlap_symmetric_paired.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        using namespace chem::kumo;
        using Clock = std::chrono::steady_clock;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd total = molecule.density_matrix(Spin::Total);
            const Eigen::MatrixXd spin = molecule.density_matrix(Spin::Difference);
            const Eigen::MatrixXd alpha = 0.5 * (total + spin);
            const Eigen::MatrixXd beta = 0.5 * (total - spin);
            const std::vector<Eigen::MatrixXd> densities{alpha, beta, total, spin};
            const auto t0 = Clock::now();
            std::vector<Eigen::MatrixXd> references;
            for (const auto& p : densities)
                references.push_back(local_density_overlap_symmetric_paired(basis, p));
            const double separate_seconds = std::chrono::duration<double>(Clock::now() - t0).count();
            const auto multi = local_density_overlap_multichannel(basis, densities);
            if (multi.matrices.size() != densities.size())
                throw std::runtime_error("multichannel result count mismatch");
            double max_error = 0.0;
            for (std::size_t k = 0; k < densities.size(); ++k) {
                const double scale = std::max(1.0, references[k].cwiseAbs().maxCoeff());
                const double error = (multi.matrices[k] - references[k]).cwiseAbs().maxCoeff();
                const double symmetry = (multi.matrices[k] - multi.matrices[k].transpose()).cwiseAbs().maxCoeff();
                if (!multi.matrices[k].allFinite() || error > 1e-11 * scale || symmetry > 1e-12)
                    throw std::runtime_error("multichannel reference mismatch");
                max_error = std::max(max_error, error);
            }
            const double linearity = (multi.matrices[0] + multi.matrices[1] - multi.matrices[2]).cwiseAbs().maxCoeff();
            const double spin_linearity = (multi.matrices[0] - multi.matrices[1] - multi.matrices[3]).cwiseAbs().maxCoeff();
            const double scale = std::max(1.0, multi.matrices[2].cwiseAbs().maxCoeff());
            if (linearity > 1e-11 * scale || spin_linearity > 1e-11 * scale)
                throw std::runtime_error("multichannel linearity mismatch");
            std::cout << "system=" << argv[arg]
                      << " channels=" << densities.size()
                      << " shell_pairs=" << multi.statistics.shell_pairs
                      << " quartets=" << multi.statistics.evaluated_quartets
                      << " max_abs=" << max_error
                      << " linearity=" << linearity
                      << " spin_linearity=" << spin_linearity
                      << " separate_seconds=" << separate_seconds
                      << " multichannel_seconds=" << multi.statistics.total_seconds
                      << " quartet_seconds=" << multi.statistics.quartet_seconds
                      << " accumulation_seconds=" << multi.statistics.accumulation_seconds
                      << " speedup=" << separate_seconds / multi.statistics.total_seconds << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
