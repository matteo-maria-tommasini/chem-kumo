#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
struct Stats {
    std::size_t pairs = 0;
    std::size_t components = 0;
    long double primitive_terms = 0;
    double elapsed = 0;
    double largest = 0;
};
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-local-ao-norm-profile file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(9);
        for (int arg = 1; arg < argc; ++arg) {
            const std::string path = argv[arg];
            const BasisSet basis{FormattedCheckpoint(path)};
            const auto shells = local_shell_components(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            for (std::size_t a = 0; a < shells.size(); ++a)
                offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
            Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(offsets.back(), offsets.back());
            std::map<std::string, Stats> groups;
            std::size_t count = 0;
            const auto all_start = Clock::now();
            for (std::size_t a = 0; a < shells.size(); ++a) {
                for (std::size_t b = a; b < shells.size(); ++b) {
                    const auto& sa = shells[a];
                    const auto& sb = shells[b];
                    const bool same_center = sa[0][0].center == sb[0][0].center;
                    const std::string category = std::to_string(sa.size()) + "x" +
                        std::to_string(sb.size()) + (same_center ? ":same" : ":different");
                    auto& st = groups[category];
                    long double terms = 0;
                    for (const auto& ca : sa)
                        for (const auto& cb : sb)
                            terms += static_cast<long double>(ca.size()) * ca.size() *
                                     static_cast<long double>(cb.size()) * cb.size();
                    const auto start = Clock::now();
                    const auto q = local_shell_quartet({&sa, &sb, &sa, &sb});
                    for (std::size_t i = 0; i < sa.size(); ++i)
                        for (std::size_t j = 0; j < sb.size(); ++j) {
                            const double squared = q.at(i, j, i, j);
                            if (!std::isfinite(squared) || squared < -1e-12 * std::max(1.0, std::abs(squared)))
                                throw std::runtime_error("invalid AO product norm squared");
                            const double v = std::sqrt(std::max(0.0, squared));
                            const Eigen::Index mu = offsets[a] + static_cast<Eigen::Index>(i);
                            const Eigen::Index nu = offsets[b] + static_cast<Eigen::Index>(j);
                            norms(mu, nu) = norms(nu, mu) = v;
                        }
                    const double dt = seconds(start);
                    ++count;
                    ++st.pairs;
                    st.components += sa.size() * sb.size();
                    st.primitive_terms += terms;
                    st.elapsed += dt;
                    st.largest = std::max(st.largest, dt);
                }
            }
            const double total = seconds(all_start);
            std::cout << "system=" << path << " nao=" << norms.rows()
                      << " shells=" << shells.size() << " shell_pairs=" << count
                      << " profile_seconds=" << total
                      << " norm_checksum=" << norms.sum() << '\n';
            for (const auto& [category, st] : groups) {
                std::cout << " category=" << category
                          << " pairs=" << st.pairs
                          << " ao_components=" << st.components
                          << " primitive_terms=" << static_cast<double>(st.primitive_terms)
                          << " seconds=" << st.elapsed
                          << " time_fraction=" << st.elapsed / total
                          << " mean_us=" << st.elapsed * 1e6 / st.pairs
                          << " max_pair_ms=" << st.largest * 1e3 << '\n';
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "AO norm profile: " << e.what() << '\n';
        return 1;
    }
}
