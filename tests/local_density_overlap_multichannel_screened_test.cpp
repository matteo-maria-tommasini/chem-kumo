#include "chem/kumo/local_density_overlap_multichannel_screened.hpp"
#include "chem/kumo/local_density_overlap_streaming.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd total = molecule.density_matrix(Spin::Total);
            const Eigen::MatrixXd spin = molecule.density_matrix(Spin::Difference);
            const Eigen::MatrixXd alpha = 0.5 * (total + spin);
            const Eigen::MatrixXd beta = 0.5 * (total - spin);
            const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(total.rows(), total.cols());
            const std::vector<Eigen::MatrixXd> densities{alpha, beta, spin, zero};
            const auto exact = local_density_overlap_multichannel(basis, densities);
            for (double factor : {0.0, 1.0}) {
                const std::vector<double> eps{factor*1e-8, factor*1e-10, factor*1e-6, 0.0};
                const auto result = local_density_overlap_multichannel_screened(basis, densities, eps, 1024*1024);
                if (result.matrices.size() != densities.size() ||
                    result.statistics.evaluated_quartets + result.statistics.skipped_quartets !=
                        result.statistics.total_quartets)
                    throw std::runtime_error("screened multichannel: quartet accounting mismatch");
                for (std::size_t k = 0; k < densities.size(); ++k) {
                    const double scale = std::max(1.0, exact.matrices[k].cwiseAbs().maxCoeff());
                    const double err = (result.matrices[k] - exact.matrices[k]).cwiseAbs().maxCoeff();
                    const double bound = result.statistics.max_error_bounds[k];
                    if (!result.matrices[k].allFinite() ||
                        err > eps[k] + 2e-11*scale || bound > eps[k] + 1e-18 ||
                        (result.matrices[k] - result.matrices[k].transpose()).cwiseAbs().maxCoeff() > 1e-11*scale)
                        throw std::runtime_error("screened multichannel: accuracy or budget failure");
                    std::cout << " channel=" << k << " epsilon=" << eps[k]
                              << " max_abs=" << err << " max_bound=" << bound;
                }
                std::cout << " system=" << argv[arg]
                          << " evaluated=" << result.statistics.evaluated_quartets
                          << " skipped=" << result.statistics.skipped_quartets
                          << " active=" << result.statistics.active_channels
                          << " peak_pair_bytes=" << result.statistics.peak_pair_bytes
                          << " total_seconds=" << result.statistics.total_seconds << '\n';
            }
            const auto zeros = local_density_overlap_multichannel_screened(basis, {zero, zero}, {0.0, 1e-8}, 0);
            if (zeros.statistics.evaluated_quartets || zeros.statistics.active_channels ||
                zeros.statistics.skipped_quartets != zeros.statistics.total_quartets)
                throw std::runtime_error("screened multichannel: all-zero fast path failed");
            bool rejected = false;
            try { (void)local_density_overlap_multichannel_screened(basis, {alpha}, {}, 1024*1024); }
            catch (const std::invalid_argument&) { rejected = true; }
            if (!rejected) throw std::runtime_error("screened multichannel: budget mismatch accepted");
            rejected = false;
            try { (void)local_density_overlap_multichannel_screened(basis, {alpha}, {-1.0}, 1024*1024); }
            catch (const std::invalid_argument&) { rejected = true; }
            if (!rejected) throw std::runtime_error("screened multichannel: negative budget accepted");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
