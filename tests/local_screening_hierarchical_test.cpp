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
#include <stdexcept>
#include <vector>

namespace {
struct Pair {
    std::size_t a, b;
    long double bmax = 0.0L;
    long double dsum = 0.0L;
    long double work = 0.0L;
};
struct Result {
    std::size_t skipped = 0;
    long double max_bound = 0.0L;
    long double saved_work = 0.0L;
    long double total_work = 0.0L;
    double seconds = 0.0;
};

// Reference streaming selector (same canonical traversal as Stage 03C.11d).
// Deterministic, strictly streaming: O(number of shell pairs) auxiliary memory.
Result select_streaming(const std::vector<Pair>& pairs, long double epsilon) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    const std::size_t n = pairs.size();
    std::vector<long double> spent(n, 0.0L);
    Result out;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& x = pairs[i];
        const long double mi = x.a == x.b ? 1.0L : 2.0L;
        for (std::size_t j = i; j < n; ++j) {
            const auto& y = pairs[j];
            const long double mj = y.a == y.b ? 1.0L : 2.0L;
            const long double to_i = mj * x.bmax * y.dsum;
            const long double to_j = i == j ? 0.0L : mi * y.bmax * x.dsum;
            const long double work = x.work * y.work;
            out.total_work += work;
            if (spent[i] + to_i > epsilon ||
                (i != j && spent[j] + to_j > epsilon)) continue;
            spent[i] += to_i;
            if (i != j) spent[j] += to_j;
            ++out.skipped;
            out.saved_work += work;
        }
    }
    if (!spent.empty()) out.max_bound = *std::max_element(spent.begin(), spent.end());
    out.seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}

// Hierarchical blocks of consecutive source shell pairs. An entire block is
// skipped only when its *aggregate* bound fits the destination budget and
// every reciprocal destination budget. Otherwise descend to individual pairs.
// No quartets or skip masks are stored. This is a conservative block hierarchy,
// not yet a spatial tree.
struct HierResult : Result {
    std::size_t block_skips = 0;
    std::size_t block_accepted_pairs = 0;
};
HierResult select_hierarchical(const std::vector<Pair>& pairs,
                               long double epsilon, std::size_t width) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    if (width == 0) throw std::runtime_error("block width must be positive");
    const std::size_t n = pairs.size();
    std::vector<long double> spent(n, 0.0L);
    HierResult out;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& x = pairs[i];
        const long double mi = x.a == x.b ? 1.0L : 2.0L;
        for (std::size_t begin = i; begin < n;) {
            const std::size_t end = begin + std::min(width, n - begin);
            long double to_i_sum = 0.0L;
            bool reciprocal_ok = true;
            long double block_work = 0.0L;
            for (std::size_t j = begin; j < end; ++j) {
                const auto& y = pairs[j];
                const long double mj = y.a == y.b ? 1.0L : 2.0L;
                to_i_sum += mj * x.bmax * y.dsum;
                if (j != i && spent[j] + mi*y.bmax*x.dsum > epsilon)
                    reciprocal_ok = false;
                block_work += x.work * y.work;
            }
            if (reciprocal_ok && spent[i] + to_i_sum <= epsilon) {
                spent[i] += to_i_sum;
                for (std::size_t j = begin; j < end; ++j) {
                    if (j != i) spent[j] += mi*pairs[j].bmax*x.dsum;
                }
                out.skipped += end - begin;
                out.saved_work += block_work;
                out.total_work += block_work;
                ++out.block_skips;
                out.block_accepted_pairs += end - begin;
            } else {
                for (std::size_t j = begin; j < end; ++j) {
                    const auto& y = pairs[j];
                    const long double mj = y.a == y.b ? 1.0L : 2.0L;
                    const long double to_i = mj*x.bmax*y.dsum;
                    const long double to_j = j == i ? 0.0L : mi*y.bmax*x.dsum;
                    const long double work = x.work*y.work;
                    out.total_work += work;
                    if (spent[i] + to_i > epsilon ||
                        (j != i && spent[j] + to_j > epsilon)) continue;
                    spent[i] += to_i;
                    if (j != i) spent[j] += to_j;
                    ++out.skipped;
                    out.saved_work += work;
                }
            }
            begin = end;
        }
    }
    if (!spent.empty()) out.max_bound = *std::max_element(spent.begin(), spent.end());
    out.seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-local-screening-hierarchical-test <fchk> [<fchk> ...]\n";
        return 2;
    }
    try {
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = 1; arg < argc; ++arg) {
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const Eigen::MatrixXd norms = local_ao_product_norms(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            std::vector<long double> primitives(shells.size(), 0.0L);
            for (std::size_t a = 0; a < shells.size(); ++a) {
                offsets[a+1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
                for (const auto& component : shells[a]) primitives[a] += component.size();
            }
            for (Spin spin : {Spin::Total, Spin::Difference}) {
                const Eigen::MatrixXd density = molecule.density_matrix(spin);
                if (density.rows() != norms.rows() || density.cols() != norms.cols() ||
                    !density.allFinite()) throw std::runtime_error("invalid density matrix");
                if ((density-density.transpose()).cwiseAbs().maxCoeff() >
                    1e-10 * std::max(1.0, density.cwiseAbs().maxCoeff()))
                    throw std::runtime_error("density matrix not symmetric");
                std::vector<Pair> pairs;
                pairs.reserve(shells.size() * (shells.size() + 1) / 2);
                for (std::size_t a = 0; a < shells.size(); ++a)
                    for (std::size_t b = a; b < shells.size(); ++b) {
                        Pair p{a, b};
                        p.work = primitives[a] * primitives[b];
                        for (Eigen::Index mu = offsets[a]; mu < offsets[a+1]; ++mu)
                            for (Eigen::Index nu = offsets[b]; nu < offsets[b+1]; ++nu) {
                                const long double v = norms(mu, nu);
                                p.bmax = std::max(p.bmax, v);
                                p.dsum += std::abs(static_cast<long double>(density(mu,nu))) * v;
                            }
                        pairs.push_back(p);
                    }
                const std::size_t n = pairs.size();
                std::cout << "system=" << argv[arg]
                          << " spin=" << (spin == Spin::Total ? "total" : "difference")
                          << " nAO=" << norms.rows() << " shell_pairs=" << n
                          << " quartets=" << n*(n+1)/2
                          << " selection_aux_bytes=" << n*sizeof(long double) << '\n';
                for (long double epsilon : {1e-6L, 1e-8L, 1e-10L, 1e-12L}) {
                    const Result ref = select_streaming(pairs, epsilon);
                    for (std::size_t width : {std::size_t{8}, std::size_t{32}, std::size_t{128}}) {
                        const HierResult r = select_hierarchical(pairs, epsilon, width);
                        if (r.max_bound > epsilon || r.skipped > n*(n+1)/2)
                            throw std::runtime_error("hierarchical budget/count exceeded");
                        const long double fraction = r.total_work > 0 ? r.saved_work/r.total_work : 0;
                        std::cout << " epsilon=" << static_cast<double>(epsilon)
                                  << " block_width=" << width
                                  << " reference_skipped=" << ref.skipped
                                  << " hierarchical_skipped=" << r.skipped
                                  << " block_skips=" << r.block_skips
                                  << " block_accepted_pairs=" << r.block_accepted_pairs
                                  << " max_error_bound=" << static_cast<double>(r.max_bound)
                                  << " estimated_saved_work_fraction=" << static_cast<double>(fraction)
                                  << " selection_seconds=" << r.seconds << '\n';
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "hierarchical screening diagnostic: " << e.what() << '\n';
        return 1;
    }
}
