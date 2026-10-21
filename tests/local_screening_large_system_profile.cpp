#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <array>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
struct PairProfile {
    double importance;
    double bmax;
    double weighted_density;
    std::size_t a, b;
    double distance_bohr;
};
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-local-screening-large-system-profile [--top N] file.fchk[.gz] ...\n";
        return 2;
    }
    try {
        std::size_t top_n = 0;
        int first_file = 1;
        if (argc >= 3 && std::string(argv[1]) == "--top") {
            const std::string count_text = argv[2];
            if (count_text.empty() || count_text.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument("--top expects a nonnegative integer");
            const auto count = std::stoull(count_text);
            if (count > 10000) throw std::invalid_argument("--top must be <= 10000");
            top_n = static_cast<std::size_t>(count);
            first_file = 3;
        }
        if (first_file >= argc) throw std::invalid_argument("missing FCHK filename");
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = first_file; arg < argc; ++arg) {
            const std::string path = argv[arg];
            const auto start = Clock::now();
            const auto molecule = Molecule::from_fchk(path);
            const BasisSet basis{FormattedCheckpoint(path)};
            const auto shells = local_shell_components(basis);
            const double load_seconds = elapsed(start);
            const auto norms_start = Clock::now();
            const Eigen::MatrixXd norms = local_ao_product_norms(basis);
            const double norms_seconds = elapsed(norms_start);
            const auto density = molecule.density_matrix(Spin::Total);
            if (density.rows() != norms.rows() || density.cols() != norms.cols() ||
                !density.allFinite() || !norms.allFinite())
                throw std::runtime_error("density/norm shape or finiteness mismatch");
            const auto profile_start = Clock::now();
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            for (std::size_t a = 0; a < shells.size(); ++a)
                offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
            std::vector<std::array<double, 3>> origins;
            origins.reserve(shells.size());
            for (const auto& shell : shells) {
                if (shell.empty() || shell.front().empty())
                    throw std::runtime_error("empty shell component");
                origins.push_back(shell.front().front().center);
            }
            std::vector<PairProfile> pairs;
            pairs.reserve(shells.size() * (shells.size() + 1) / 2);
            for (std::size_t a = 0; a < shells.size(); ++a) {
                for (std::size_t b = a; b < shells.size(); ++b) {
                    double bmax = 0.0;
                    long double dsum = 0.0L;
                    for (Eigen::Index mu = offsets[a]; mu < offsets[a + 1]; ++mu)
                        for (Eigen::Index nu = offsets[b]; nu < offsets[b + 1]; ++nu) {
                            const double bound = norms(mu, nu);
                            bmax = std::max(bmax, bound);
                            dsum += std::abs(density(mu, nu)) * bound;
                        }
                    const long double mult = a == b ? 1.0L : 2.0L;
                    const double importance = static_cast<double>(mult * bmax * dsum);
                    if (!std::isfinite(importance))
                        throw std::runtime_error("nonfinite importance");
                    double distance2 = 0.0;
                    for (int axis = 0; axis < 3; ++axis) {
                        const double delta = origins[a][axis] - origins[b][axis];
                        distance2 += delta * delta;
                    }
                    pairs.push_back({importance, bmax, static_cast<double>(dsum),
                                     a, b, std::sqrt(distance2)});
                }
            }
            const double profile_seconds = elapsed(profile_start);
            std::vector<double> importance;
            importance.reserve(pairs.size());
            for (const auto& p : pairs) importance.push_back(p.importance);
            std::sort(importance.begin(), importance.end(), std::greater<double>{});
            const long double total = std::accumulate(importance.begin(), importance.end(), 0.0L);
            std::size_t positive = 0;
            for (double v : importance) if (v > 0.0) ++positive;
            std::cout << "system=" << path
                      << " nao=" << norms.rows()
                      << " shells=" << shells.size()
                      << " shell_pairs=" << pairs.size()
                      << " formal_quartets=" << (static_cast<unsigned long long>(pairs.size()) * (pairs.size() + 1ULL) / 2ULL)
                      << " positive_importances=" << positive
                      << " norms_mib=" << (static_cast<double>(norms.size()) * sizeof(double) / 1048576.0)
                      << " load_seconds=" << load_seconds
                      << " norms_seconds=" << norms_seconds
                      << " profile_seconds=" << profile_seconds << '\n';
            for (double fraction : {0.001, 0.01, 0.05, 0.10, 0.50}) {
                const auto count = std::min(importance.size(), std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(fraction * importance.size()))));
                const long double selected = std::accumulate(importance.begin(), importance.begin() + static_cast<std::ptrdiff_t>(count), 0.0L);
                std::cout << " top_fraction=" << fraction << " top_pairs=" << count
                          << " importance_mass_fraction=" << (total > 0.0L ? static_cast<double>(selected / total) : 0.0)
                          << '\n';
            }
            if (top_n != 0) {
                std::vector<std::size_t> order(pairs.size());
                std::iota(order.begin(), order.end(), std::size_t{0});
                const auto compare = [&](std::size_t x, std::size_t y) {
                    if (pairs[x].importance != pairs[y].importance)
                        return pairs[x].importance > pairs[y].importance;
                    return x < y;
                };
                const auto count = std::min(top_n, order.size());
                std::partial_sort(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(count),
                                  order.end(), compare);
                long double cumulative = 0.0L;
                for (std::size_t rank = 0; rank < count; ++rank) {
                    const auto& pair = pairs[order[rank]];
                    cumulative += pair.importance;
                    std::cout << " dominant_rank=" << rank + 1
                              << " shell_a=" << pair.a << " shell_b=" << pair.b
                              << " center_distance_bohr=" << pair.distance_bohr
                              << " components_a=" << shells[pair.a].size()
                              << " components_b=" << shells[pair.b].size()
                              << " bmax=" << pair.bmax
                              << " weighted_density=" << pair.weighted_density
                              << " importance=" << pair.importance
                              << " cumulative_mass_fraction="
                              << (total > 0.0L ? static_cast<double>(cumulative / total) : 0.0)
                              << '\n';
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "large-system profile: " << e.what() << '\n';
        return 1;
    }
}
