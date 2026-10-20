#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/basis.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const chem::kumo::FormattedCheckpoint fchk(argv[arg]);
            const chem::kumo::BasisSet basis(fchk);
            const auto shells = chem::kumo::local_shell_components(basis);
            const auto norms = chem::kumo::local_ao_product_norms(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            for (std::size_t a = 0; a < shells.size(); ++a)
                offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
            double max_violation = 0.0;
            std::size_t checked = 0;
            const std::size_t limit = std::min<std::size_t>(shells.size(), 6);
            for (std::size_t a = 0; a < limit; ++a)
                for (std::size_t b = 0; b < limit; ++b)
                    for (std::size_t c = 0; c < limit; ++c)
                        for (std::size_t d = 0; d < limit; ++d) {
                            const auto q = chem::kumo::local_shell_quartet({&shells[a], &shells[b], &shells[c], &shells[d]});
                            for (std::size_t i = 0; i < q.dimensions[0]; ++i)
                                for (std::size_t j = 0; j < q.dimensions[1]; ++j)
                                    for (std::size_t k = 0; k < q.dimensions[2]; ++k)
                                        for (std::size_t l = 0; l < q.dimensions[3]; ++l) {
                                            const double bound = norms(offsets[a]+i, offsets[b]+j) *
                                                norms(offsets[c]+k, offsets[d]+l);
                                            const double value = std::abs(q.at(i,j,k,l));
                                            const double violation = value - bound;
                                            max_violation = std::max(max_violation, violation);
                                            if (violation > 2e-11 * std::max({1.0, value, bound}))
                                                throw std::runtime_error("Cauchy-Schwarz violation");
                                            ++checked;
                                        }
                        }
            std::cout << std::setprecision(12) << argv[arg]
                      << " nAO=" << norms.rows() << " quartets=" << checked
                      << " max_positive_violation=" << max_violation << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
