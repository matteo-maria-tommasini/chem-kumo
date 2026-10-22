#include "chem/kumo/local_w_dispatch.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        using namespace chem::kumo;
        std::size_t checks = 0;
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const auto total = molecule.density_matrix(Spin::Total);
            const auto spin = molecule.density_matrix(Spin::Difference);
            const std::vector<Eigen::MatrixXd> densities{total, spin};
            const auto reference = local_density_overlap_multichannel(basis, densities);
            const auto via_reference = local_w_evaluate(basis, densities);
            if (via_reference.backend != LocalWBackend::Reference ||
                via_reference.matrices.size() != densities.size())
                throw std::runtime_error("W reference dispatch metadata mismatch");
            for (std::size_t k = 0; k < densities.size(); ++k) {
                if (!via_reference.matrices[k].isApprox(reference.matrices[k], 0.0))
                    throw std::runtime_error("W reference dispatch mismatch");
                ++checks;
            }
            for (const auto backend : {LocalWBackend::Screened, LocalWBackend::Hierarchical}) {
                LocalWOptions options;
                options.backend = backend;
                options.error_budgets = {1e-8, 1e-10};
                options.screening_order = MultichannelScreeningOrder::MaximumImportance;
                const auto result = local_w_evaluate(basis, densities, options);
                if (result.backend != backend ||
                    result.statistics.evaluated_quartets + result.statistics.skipped_quartets !=
                        result.statistics.total_quartets)
                    throw std::runtime_error("W dispatch quartet accounting mismatch");
                for (std::size_t k = 0; k < densities.size(); ++k) {
                    const double scale = std::max(1.0, reference.matrices[k].cwiseAbs().maxCoeff());
                    const double error = (result.matrices[k] - reference.matrices[k]).cwiseAbs().maxCoeff();
                    if (!result.matrices[k].allFinite() ||
                        error > options.error_budgets[k] + 2e-11 * scale ||
                        result.statistics.max_error_bounds[k] > options.error_budgets[k] + 1e-18)
                        throw std::runtime_error("W dispatch error budget mismatch");
                    ++checks;
                }
            }
            LocalWOptions bad;
            bad.backend = LocalWBackend::Screened;
            try {
                (void)local_w_evaluate(basis, densities, bad);
                throw std::runtime_error("W dispatch accepted missing budgets");
            } catch (const std::invalid_argument&) { ++checks; }
            std::cout << "system=" << argv[arg] << " status=OK\n";
        }
        std::cout << "checks=" << checks << " status=OK\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
