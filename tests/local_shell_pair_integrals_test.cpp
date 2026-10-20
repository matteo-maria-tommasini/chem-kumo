#include "chem/kumo/local_shell_pair_integrals.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    using namespace chem::kumo;
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const auto limit = std::min<std::size_t>(shells.size(), 8);
            using Clock = std::chrono::steady_clock;
            const auto prep_begin = Clock::now();
            std::vector<LocalPreparedShellPair> pairs;
            pairs.reserve(limit * limit);
            std::size_t prepared_primitives = 0;
            for (std::size_t a = 0; a < limit; ++a)
                for (std::size_t b = 0; b < limit; ++b) {
                    pairs.push_back(prepare_local_shell_pair(shells[a], shells[b]));
                    prepared_primitives += pairs.back().pair_count();
                }
            const auto prep_end = Clock::now();
            double max_abs = 0.0, max_scaled = 0.0;
            std::size_t quartets = 0, elements = 0;
            double direct_seconds = 0.0, paired_seconds = 0.0;
            volatile double checksum = 0.0;
            // A fixed, diverse subset of shell quartets; full AO-component
            // comparison for each chosen quartet (including pure d components).
            for (std::size_t a = 0; a < limit; ++a)
                for (std::size_t b = 0; b < limit; ++b)
                    for (std::size_t c = 0; c < limit; ++c) {
                        const std::size_t d = (a + 2*b + 3*c) % limit;
                        const auto start_direct = Clock::now();
                        const auto reference = local_shell_quartet(
                            {&shells[a], &shells[b], &shells[c], &shells[d]});
                        const auto end_direct = Clock::now();
                        const auto actual = local_shell_quartet_paired(
                            pairs[a*limit+b], pairs[c*limit+d]);
                        const auto end_pair = Clock::now();
                        direct_seconds += std::chrono::duration<double>(end_direct-start_direct).count();
                        paired_seconds += std::chrono::duration<double>(end_pair-end_direct).count();
                        if (reference.dimensions != actual.dimensions ||
                            reference.values.size() != actual.values.size())
                            throw std::runtime_error("paired quartet shape mismatch");
                        for (std::size_t n = 0; n < reference.values.size(); ++n) {
                            const double x = reference.values[n], y = actual.values[n];
                            const double err = std::abs(x-y);
                            const double scale = std::max({1.0,std::abs(x),std::abs(y)});
                            max_abs = std::max(max_abs, err);
                            max_scaled = std::max(max_scaled, err/scale);
                            if (!std::isfinite(x) || !std::isfinite(y) || err > 1e-12*scale)
                                throw std::runtime_error("paired shell quartet mismatch");
                            checksum = checksum + y;
                            ++elements;
                        }
                        ++quartets;
                    }
            std::cout << std::setprecision(12)
                      << "system=" << argv[arg] << " shells_tested=" << limit
                      << " quartets=" << quartets << " elements=" << elements
                      << " prepared_primitive_pairs=" << prepared_primitives
                      << " pair_storage_bytes_lower_bound="
                      << prepared_primitives*sizeof(GaussianPairIntermediate)
                      << " preparation_seconds="
                      << std::chrono::duration<double>(prep_end-prep_begin).count()
                      << " direct_seconds=" << direct_seconds
                      << " paired_seconds=" << paired_seconds
                      << " speedup_kernel=" << direct_seconds/paired_seconds
                      << " max_abs=" << max_abs << " max_scaled=" << max_scaled
                      << " checksum=" << checksum << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
