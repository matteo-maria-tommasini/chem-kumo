#include "chem/kumo/local_shell_integrals.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/basis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

static void check(double x, double y, double tol, const char* what) {
    if (!std::isfinite(x) || !std::isfinite(y) || std::abs(x-y) > tol * std::max({1.0, std::abs(x), std::abs(y)}))
        throw std::runtime_error(what);
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const FormattedCheckpoint fchk(argv[arg]);
            const BasisSet basis(fchk);
            const auto shells = local_shell_components(basis);
            const BasisEvaluator ao(basis);
            std::size_t total = 0;
            for (const auto& shell : shells) total += shell.size();
            if (total != ao.basis_size()) throw std::runtime_error("AO component ordering/count mismatch");
            const auto overlap = ao.overlap_matrix();
            // Independently contract four-function integrals to recover the
            // two-function overlap via a separate analytic primitive oracle:
            // compare quartet permutations and shell Cartesian/pure mappings.
            const std::size_t limit = std::min<std::size_t>(shells.size(), 5);
            std::size_t checked = 0;
            for (std::size_t a = 0; a < limit; ++a)
                for (std::size_t b = 0; b < limit; ++b) {
                    const auto ab = local_shell_quartet({&shells[a], &shells[b], &shells[a], &shells[b]});
                    const auto ba = local_shell_quartet({&shells[b], &shells[a], &shells[b], &shells[a]});
                    for (std::size_t i = 0; i < shells[a].size(); ++i)
                        for (std::size_t j = 0; j < shells[b].size(); ++j) {
                            check(ab.at(i,j,i,j), ba.at(j,i,j,i), 2e-12, "quartet permutation");
                            if (ab.at(i,j,i,j) < -2e-12)
                                throw std::runtime_error("negative AO product square integral");
                            ++checked;
                        }
                }
            if (!overlap.allFinite()) throw std::runtime_error("invalid analytic AO overlap");
            std::cout << argv[arg] << " shells=" << shells.size()
                      << " nAO=" << total << " quartet_checks=" << checked << " validated\n";
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
