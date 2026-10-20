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
#include <stdexcept>
#include <vector>
#include <utility>

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

// Reorder shell pairs by ascending importance s_i = m_i b_i d_i.
// The blocked variant traverses the same sorted order in contiguous chunks.
// It is deliberately a control benchmark: no unsupported claim of block pruning.
Result select_ordered(const std::vector<Pair>& pairs, long double epsilon,
                      std::size_t block_width) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    const std::size_t n = pairs.size();
    std::vector<std::size_t> order(n);
    for (std::size_t i=0; i<n; ++i) order[i]=i;
    const auto importance = [&](std::size_t i) {
        const auto& p=pairs[i];
        return (p.a==p.b ? 1.0L : 2.0L)*p.bmax*p.dsum;
    };
    std::sort(order.begin(),order.end(),[&](std::size_t i,std::size_t j) {
        const auto si=importance(i), sj=importance(j);
        return si<sj || (si==sj && i<j);
    });
    std::vector<long double> spent(n,0.0L);
    Result out;
    // The blocks only change loop organization, not the visit order.
    for (std::size_t ib=0; ib<n; ib+=block_width) {
        const std::size_t ie=std::min(n,ib+block_width);
        for (std::size_t ii=ib; ii<ie; ++ii) {
            const std::size_t i=order[ii];
            const auto& x=pairs[i];
            const long double mi=x.a==x.b ? 1.0L:2.0L;
            for (std::size_t jj=ii; jj<n; ++jj) {
                const std::size_t j=order[jj];
                const auto& y=pairs[j];
                const long double mj=y.a==y.b ? 1.0L:2.0L;
                const long double ei=mj*x.bmax*y.dsum;
                const long double ej=i==j ? 0.0L:mi*y.bmax*x.dsum;
                const long double work=x.work*y.work;
                out.total_work+=work;
                if (spent[i]+ei>epsilon || (i!=j && spent[j]+ej>epsilon)) continue;
                spent[i]+=ei;
                if (i!=j) spent[j]+=ej;
                ++out.skipped;
                out.saved_work+=work;
            }
        }
    }
    if (!spent.empty()) out.max_bound=*std::max_element(spent.begin(),spent.end());
    out.seconds=std::chrono::duration<double>(Clock::now()-start).count();
    return out;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-local-screening-lazy-test <fchk> [<fchk> ...]\n";
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
                          << " reference_aux_bytes=" << n*sizeof(long double) << '\n';
                for (long double epsilon : {1e-6L, 1e-8L, 1e-10L, 1e-12L}) {
                    const Result ref = select_streaming(pairs, epsilon);
                    const Result imp = select_ordered(pairs, epsilon, 1);
                    const Result blk = select_ordered(pairs, epsilon, 32);
                    if (ref.max_bound > epsilon + 1e-17L ||
                        imp.max_bound > epsilon + 1e-17L ||
                        blk.max_bound > epsilon + 1e-17L)
                        throw std::runtime_error("error budget exceeded");
                    for (const auto& [name, r] :
                         std::array<std::pair<const char*,Result>,3>{{
                           {"canonical",ref}, {"importance",imp}, {"importance_blocked",blk}}}) {
                        const long double fraction = r.total_work > 0 ?
                            r.saved_work / r.total_work : 0;
                        std::cout << " epsilon=" << static_cast<double>(epsilon)
                                  << " strategy=" << name
                                  << " skipped=" << r.skipped
                                  << " max_error_bound=" << static_cast<double>(r.max_bound)
                                  << " saved_work_fraction=" << static_cast<double>(fraction)
                                  << " selection_seconds=" << r.seconds << '\n';
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "importance screening diagnostic: " << e.what() << '\n';
        return 1;
    }
}
