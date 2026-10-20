#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
struct Pair {
    std::size_t a, b;
    long double bmax = 0.0L;
    long double work = 0.0L;
    std::vector<long double> dsum;
};
struct Result {
    std::size_t skipped = 0;
    long double saved_work = 0.0L;
    long double total_work = 0.0L;
    std::vector<long double> bounds;
    double ordering_seconds = 0.0;
    double selection_seconds = 0.0;
};
enum class Order { Canonical, Maximum, Sum };

Result select(const std::vector<Pair>& pairs,
              const std::vector<long double>& budgets, Order mode) {
    using Clock = std::chrono::steady_clock;
    const std::size_t n = pairs.size(), channels = budgets.size();
    const auto begin = Clock::now();
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (mode != Order::Canonical) {
        std::vector<long double> scores(n, 0.0L);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& p = pairs[i];
            const long double mult = p.a == p.b ? 1.0L : 2.0L;
            for (std::size_t k = 0; k < channels; ++k) {
                // A zero budget is a hard constraint, not a denominator.
                // Positive importance at zero budget is ranked last.
                const long double numerator = mult * p.bmax * p.dsum[k];
                const long double term = budgets[k] == 0.0L
                    ? (numerator == 0.0L ? 0.0L : std::numeric_limits<long double>::infinity())
                    : numerator / budgets[k];
                if (mode == Order::Maximum) scores[i] = std::max(scores[i], term);
                else scores[i] += term;
            }
        }
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return scores[a] < scores[b] || (scores[a] == scores[b] && a < b);
        });
    }
    const auto ordered = Clock::now();
    std::vector<std::vector<long double>> spent(channels, std::vector<long double>(n, 0.0L));
    Result out;
    out.bounds.assign(channels, 0.0L);
    out.ordering_seconds = std::chrono::duration<double>(ordered - begin).count();
    for (std::size_t ii = 0; ii < n; ++ii) {
        const std::size_t i = order[ii];
        const auto& x = pairs[i];
        const long double mi = x.a == x.b ? 1.0L : 2.0L;
        for (std::size_t jj = ii; jj < n; ++jj) {
            const std::size_t j = order[jj];
            const auto& y = pairs[j];
            const long double mj = y.a == y.b ? 1.0L : 2.0L;
            const long double work = x.work * y.work;
            out.total_work += work;
            bool skip = true;
            for (std::size_t k = 0; k < channels; ++k) {
                const long double ei = mj * x.bmax * y.dsum[k];
                const long double ej = i == j ? 0.0L : mi * y.bmax * x.dsum[k];
                if (spent[k][i] + ei > budgets[k] ||
                    (i != j && spent[k][j] + ej > budgets[k])) {
                    skip = false;
                    break;
                }
            }
            if (!skip) continue;
            for (std::size_t k = 0; k < channels; ++k) {
                spent[k][i] += mj * x.bmax * y.dsum[k];
                if (i != j) spent[k][j] += mi * y.bmax * x.dsum[k];
            }
            ++out.skipped;
            out.saved_work += work;
        }
    }
    for (std::size_t k = 0; k < channels; ++k)
        if (!spent[k].empty()) out.bounds[k] = *std::max_element(spent[k].begin(), spent[k].end());
    out.selection_seconds = std::chrono::duration<double>(Clock::now() - ordered).count();
    return out;
}
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-local-screening-multichannel-importance-test <fchk> [<fchk> ...]\n";
        return 2;
    }
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
            const auto shells = local_shell_components(basis);
            const auto norms = local_ao_product_norms(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            std::vector<long double> primitive_entries(shells.size(), 0.0L);
            for (std::size_t a = 0; a < shells.size(); ++a) {
                offsets[a+1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
                for (const auto& component : shells[a]) primitive_entries[a] += component.size();
            }
            if (offsets.back() != total.rows()) throw std::runtime_error("AO dimension mismatch");
            std::vector<std::size_t> active;
            for (std::size_t k = 0; k < densities.size(); ++k) {
                if (!densities[k].allFinite()) throw std::runtime_error("nonfinite density");
                if (!(densities[k].array() == 0.0).all()) active.push_back(k);
            }
            std::vector<Pair> pairs;
            for (std::size_t a = 0; a < shells.size(); ++a)
                for (std::size_t b = a; b < shells.size(); ++b) {
                    Pair p{a, b, 0.0L, 0.0L, {}};
                    p.work = primitive_entries[a] * primitive_entries[b];
                    p.dsum.assign(active.size(), 0.0L);
                    for (Eigen::Index mu = offsets[a]; mu < offsets[a+1]; ++mu)
                        for (Eigen::Index nu = offsets[b]; nu < offsets[b+1]; ++nu) {
                            const long double v = norms(mu, nu);
                            p.bmax = std::max(p.bmax, v);
                            for (std::size_t k = 0; k < active.size(); ++k)
                                p.dsum[k] += std::abs(static_cast<long double>(densities[active[k]](mu, nu))) * v;
                        }
                    pairs.push_back(std::move(p));
                }
            const std::size_t np = pairs.size();
            std::cout << "system=" << argv[arg] << " nao=" << total.rows()
                      << " shell_pairs=" << np << " quartets=" << np*(np+1)/2
                      << " active_channels=" << active.size() << '\n';
            for (const auto& all_budgets : std::array<std::array<long double,4>,4>{{
                    {{1e-8L,1e-10L,1e-6L,0.0L}},
                    {{1e-6L,1e-6L,1e-6L,0.0L}},
                    {{1e-10L,1e-10L,1e-10L,0.0L}},
                    {{0.0L,1e-8L,1e-6L,0.0L}} }}) {
                std::vector<long double> budgets;
                for (auto k : active) budgets.push_back(all_budgets[k]);
                for (const auto& [name, mode] : std::array<std::pair<const char*,Order>,3>{{
                         {"canonical",Order::Canonical}, {"maximum",Order::Maximum}, {"sum",Order::Sum}}}) {
                    const auto r = select(pairs, budgets, mode);
                    for (std::size_t k = 0; k < budgets.size(); ++k)
                        if (r.bounds[k] > budgets[k] + 1e-18L)
                            throw std::runtime_error("error budget exceeded");
                    std::cout << " budgets=";
                    for (std::size_t k = 0; k < budgets.size(); ++k)
                        std::cout << (k ? "," : "") << static_cast<double>(budgets[k]);
                    std::cout << " strategy=" << name
                              << " skipped=" << r.skipped
                              << " saved_work_fraction=" << static_cast<double>(r.total_work ? r.saved_work/r.total_work : 0.0L)
                              << " max_bounds=";
                    for (std::size_t k = 0; k < r.bounds.size(); ++k)
                        std::cout << (k ? "," : "") << static_cast<double>(r.bounds[k]);
                    std::cout << " ordering_seconds=" << r.ordering_seconds
                              << " selection_seconds=" << r.selection_seconds << '\n';
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "multichannel importance diagnostic: " << e.what() << '\n';
        return 1;
    }
}
