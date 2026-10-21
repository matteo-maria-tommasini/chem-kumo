#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_integrals.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/basis.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

int main(int argc, char** argv) {
    using namespace chem::kumo;
    if (argc < 2) {
        std::cerr << "usage: norm-regression file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const Eigen::MatrixXd actual = local_ao_product_norms(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            for (std::size_t a = 0; a < shells.size(); ++a)
                offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
            const bool exhaustive = shells.size() <= 32;
            std::vector<std::pair<std::size_t, std::size_t>> selected;
            if (exhaustive) {
                for (std::size_t a = 0; a < shells.size(); ++a)
                    for (std::size_t b = a; b < shells.size(); ++b)
                        selected.emplace_back(a, b);
            } else {
                // Cover same-shell, distant-shell, and all angular component sizes.
                // Reference quartets are expensive for f-f shells: sample deliberately.
                std::vector<std::size_t> representatives;
                for (std::size_t components : {1u, 3u, 5u, 7u}) {
                    for (std::size_t a = 0; a < shells.size(); ++a)
                        if (shells[a].size() == components) {
                            representatives.push_back(a);
                            break;
                        }
                }
                for (const auto a : representatives)
                    for (const auto b : representatives)
                        if (a <= b) selected.emplace_back(a, b);
                for (const auto a : representatives)
                    selected.emplace_back(a, a);
                for (std::size_t a = 0; a < shells.size(); a += std::max<std::size_t>(1, shells.size()/12)) {
                    selected.emplace_back(a, a);
                    selected.emplace_back(a, shells.size() - 1);
                }
                std::sort(selected.begin(), selected.end());
                selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
            }
            double max_abs = 0.0, max_scaled = 0.0;
            std::size_t checked = 0;
            for (const auto [a, b] : selected) {
                const auto q = local_shell_quartet({&shells[a], &shells[b], &shells[a], &shells[b]});
                for (std::size_t i = 0; i < shells[a].size(); ++i)
                    for (std::size_t j = 0; j < shells[b].size(); ++j) {
                        const double squared = q.at(i, j, i, j);
                        if (!std::isfinite(squared) || squared < -1e-12 * std::max(1.0, std::abs(squared)))
                            throw std::runtime_error("invalid reference squared norm");
                        const double reference = std::sqrt(std::max(0.0, squared));
                        const double value = actual(offsets[a] + static_cast<Eigen::Index>(i),
                                                    offsets[b] + static_cast<Eigen::Index>(j));
                        const double error = std::abs(value - reference);
                        const double scaled = error / std::max({1.0, std::abs(value), std::abs(reference)});
                        max_abs = std::max(max_abs, error);
                        max_scaled = std::max(max_scaled, scaled);
                        if (!std::isfinite(value) || scaled > 1e-11)
                            throw std::runtime_error("AO product norm regression mismatch");
                        ++checked;
                    }
            }
            const double symmetry = (actual - actual.transpose()).cwiseAbs().maxCoeff();
            if (symmetry > 1e-13) throw std::runtime_error("AO norm matrix not symmetric");
            std::cout << std::setprecision(12) << "system=" << argv[arg]
                      << " exhaustive=" << exhaustive << " shell_pairs=" << selected.size()
                      << " ao_norms_checked=" << checked << " max_abs=" << max_abs
                      << " max_scaled=" << max_scaled << " symmetry=" << symmetry << " OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
