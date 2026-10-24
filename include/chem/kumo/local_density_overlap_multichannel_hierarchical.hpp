#pragma once

#include "chem/kumo/local_density_overlap_multichannel_screened.hpp"
#include "chem/kumo/local_pair_bounds.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

struct LocalHierarchicalScreeningStatistics {
    std::size_t shell_pairs = 0;
    std::size_t total_quartets = 0;
    std::size_t evaluated_quartets = 0;
    std::size_t skipped_quartets = 0;
    std::size_t prepared_pair_instances = 0;
    std::size_t peak_pair_bytes = 0;
    std::size_t active_channels = 0;
    std::vector<double> max_error_bounds;
    std::size_t binary_nodes = 0;
    std::size_t candidate_visits = 0;
    std::size_t cutoff_skipped_quartets = 0;
    double ordering_seconds = 0.0;
    double norm_seconds = 0.0;
    double contraction_seconds = 0.0;
    double total_seconds = 0.0;
};

struct LocalHierarchicalDensityOverlap {
    std::vector<Eigen::MatrixXd> matrices;
    LocalHierarchicalScreeningStatistics statistics;
};

// Certified directed suffixes with reciprocal intersection. Only surviving
// quartets are enumerated; the memory cap applies to Gaussian pair vectors.
[[nodiscard]] inline LocalHierarchicalDensityOverlap
local_density_overlap_multichannel_hierarchical(
    const BasisSet& basis,
    const std::vector<Eigen::MatrixXd>& densities,
    const std::vector<double>& epsilons,
    std::size_t memory_limit_bytes,
    MultichannelScreeningOrder order = MultichannelScreeningOrder::MaximumImportance,
    bool screening_only = false) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    if (densities.size() != epsilons.size())
        throw std::invalid_argument("multichannel screening: budget count mismatch");
    for (double e : epsilons)
        if (!(e >= 0.0) || !std::isfinite(e))
            throw std::invalid_argument("multichannel screening: invalid budget");

    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao = offsets.back();
    std::vector<std::size_t> active;
    for (std::size_t k = 0; k < densities.size(); ++k) {
        const auto& p = densities[k];
        if (p.rows() != nao || p.cols() != nao || !p.allFinite())
            throw std::invalid_argument("multichannel screening: invalid density");
        if (nao && (p - p.transpose()).cwiseAbs().maxCoeff() >
                1e-10 * std::max(1.0, p.cwiseAbs().maxCoeff()))
            throw std::invalid_argument("multichannel screening: asymmetric density");
        if (nao && !(p.array() == 0.0).all()) active.push_back(k);
    }
    LocalHierarchicalDensityOverlap out;
    out.matrices.reserve(densities.size());
    for (std::size_t k = 0; k < densities.size(); ++k)
        out.matrices.emplace_back(Eigen::MatrixXd::Zero(nao, nao));
    auto& stats = out.statistics;
    stats.active_channels = active.size();
    stats.max_error_bounds.assign(densities.size(), 0.0);
    struct Pair { std::size_t a, b; long double bmax = 0.0L; };
    std::vector<Pair> pairs;
    for (std::size_t a = 0; a < shells.size(); ++a)
        for (std::size_t b = a; b < shells.size(); ++b)
            pairs.push_back({a, b});
    const std::size_t np = pairs.size();
    stats.shell_pairs = np;
    stats.total_quartets = np * (np + 1) / 2;
    if (active.empty()) {
        stats.skipped_quartets = stats.total_quartets;
        stats.total_seconds = std::chrono::duration<double>(Clock::now() - start).count();
        return out;
    }

    const auto norm_start = Clock::now();
    const Eigen::MatrixXd norms = local_ao_product_norms(basis);
    std::vector<std::vector<long double>> dsum(densities.size(), std::vector<long double>(np, 0.0L));
    for (std::size_t i = 0; i < np; ++i) {
        auto& pair = pairs[i];
        for (Eigen::Index mu = offsets[pair.a]; mu < offsets[pair.a+1]; ++mu)
            for (Eigen::Index nu = offsets[pair.b]; nu < offsets[pair.b+1]; ++nu) {
                const long double b = norms(mu, nu);
                pair.bmax = std::max(pair.bmax, b);
                for (std::size_t k : active)
                    dsum[k][i] += std::abs(static_cast<long double>(densities[k](mu, nu))) * b;
            }
    }
    stats.norm_seconds = std::chrono::duration<double>(Clock::now() - norm_start).count();
    const auto ordering_start = Clock::now();
    // Descending importance, as in the validated bidirectional diagnostic.
    std::vector<std::size_t> traversal(np);
    std::iota(traversal.begin(), traversal.end(), std::size_t{0});
    if (order != MultichannelScreeningOrder::MaximumImportance)
        throw std::invalid_argument("hierarchical screening: only maximum importance is supported");
    std::vector<long double> scores(np, 0.0L);
    std::vector<std::vector<long double>> weighted(active.size(),
                                                    std::vector<long double>(np, 0.0L));
    for (std::size_t i = 0; i < np; ++i) {
        const auto& p = pairs[i];
        const long double m = p.a == p.b ? 1.0L : 2.0L;
        for (std::size_t c = 0; c < active.size(); ++c) {
            const auto k = active[c];
            weighted[c][i] = m * dsum[k][i];
            const long double numerator = p.bmax * weighted[c][i];
            const long double score = epsilons[k] == 0.0
                ? (numerator == 0.0L ? 0.0L : std::numeric_limits<long double>::infinity())
                : numerator / static_cast<long double>(epsilons[k]);
            scores[i] = std::max(scores[i], score);
        }
    }
    std::sort(traversal.begin(), traversal.end(), [&](std::size_t a, std::size_t b) {
        return scores[a] > scores[b] || (scores[a] == scores[b] && a < b);
    });
    std::vector<std::vector<long double>> suffix(active.size(),
                                                  std::vector<long double>(np + 1, 0.0L));
    for (std::size_t r = np; r-- > 0;)
        for (std::size_t c = 0; c < active.size(); ++c)
            suffix[c][r] = suffix[c][r+1] + weighted[c][traversal[r]];
    std::vector<std::size_t> cut(np, np);
    for (std::size_t r = 0; r < np; ++r) {
        std::size_t lo = 0, hi = np;
        while (lo < hi) {
            const std::size_t mid = lo + (hi-lo)/2;
            bool fits = true;
            for (std::size_t c = 0; c < active.size(); ++c) {
                const auto k = active[c];
                // Conservative margin also protects the suffix roundoff.
                if (pairs[traversal[r]].bmax * suffix[c][mid] >
                    static_cast<long double>(epsilons[k]) * (1.0L - 1e-12L)) {
                    fits = false;
                    break;
                }
            }
            ++stats.binary_nodes;
            if (fits) hi = mid; else lo = mid + 1;
        }
        cut[r] = lo;
    }
    // Active set: ranks j whose reciprocal cutoff still exceeds the current i.
    // These cannot be skipped even when j is in i's certified suffix.
    std::vector<std::vector<std::size_t>> deactivate(np+1);
    for (std::size_t j = 0; j < np; ++j)
        deactivate[cut[j]].push_back(j);
    std::set<std::size_t> reciprocal_active;
    for (std::size_t j = 0; j < np; ++j)
        if (cut[j] > 0) reciprocal_active.insert(j);
    stats.ordering_seconds = std::chrono::duration<double>(Clock::now() - ordering_start).count();
    if (screening_only) {
        // Fenwick tree: count reciprocal-active ranks without enumerating quartets.
        std::vector<std::size_t> bit(np + 1, 0);
        const auto add_rank = [&](std::size_t rank, int delta) {
            for (std::size_t k = rank + 1; k <= np; k += k & (~k + 1)) {
                if (delta > 0) ++bit[k]; else --bit[k];
            }
        };
        const auto count_below = [&](std::size_t end) {
            std::size_t n = 0;
            for (std::size_t k = end; k; k -= k & (~k + 1)) n += bit[k];
            return n;
        };
        for (std::size_t j = 0; j < np; ++j)
            if (cut[j] > 0) add_rank(j, +1);
        std::size_t selected = 0;
        for (std::size_t ii = 0; ii < np; ++ii) {
            for (std::size_t j : deactivate[ii])
                if (cut[j] > 0) add_rank(j, -1);
            const std::size_t prefix_end = std::max(ii, cut[ii]);
            const std::size_t prefix = prefix_end - ii;
            const std::size_t reciprocal = count_below(np) - count_below(prefix_end);
            if (prefix > stats.total_quartets - selected ||
                reciprocal > stats.total_quartets - selected - prefix)
                throw std::overflow_error("hierarchical screening: candidate count overflow");
            selected += prefix + reciprocal;
        }
        stats.candidate_visits = selected;
        stats.evaluated_quartets = 0; // No quartet integrals were evaluated.
        stats.skipped_quartets = stats.total_quartets - selected;
        stats.cutoff_skipped_quartets = stats.skipped_quartets;
        for (std::size_t c = 0; c < active.size(); ++c) {
            const auto k = active[c];
            long double maximum = 0.0L;
            for (std::size_t r = 0; r < np; ++r)
                maximum = std::max(maximum, pairs[traversal[r]].bmax * suffix[c][cut[r]]);
            stats.max_error_bounds[k] = static_cast<double>(maximum);
        }
        stats.total_seconds = std::chrono::duration<double>(Clock::now() - start).count();
        return out;
    }

    const auto pair_bytes = [](const LocalPreparedShellPair& p) {
        std::size_t n = p.components.capacity() * sizeof(std::vector<GaussianPairIntermediate>);
        for (const auto& c : p.components) {
            if (c.capacity() > (std::numeric_limits<std::size_t>::max() - n) / sizeof(GaussianPairIntermediate))
                throw std::overflow_error("multichannel screening: pair byte overflow");
            n += c.capacity() * sizeof(GaussianPairIntermediate);
        }
        return n;
    };
    const auto estimate = [&](const Pair& p) {
        const auto& a = shells[p.a];
        const auto& b = shells[p.b];
        if (a.size() && b.size() > std::numeric_limits<std::size_t>::max() / a.size())
            throw std::overflow_error("multichannel screening: header overflow");
        const std::size_t components = a.size() * b.size();
        if (components > std::numeric_limits<std::size_t>::max() / sizeof(std::vector<GaussianPairIntermediate>))
            throw std::overflow_error("multichannel screening: header overflow");
        std::size_t bytes = components * sizeof(std::vector<GaussianPairIntermediate>);
        for (const auto& ca : a)
            for (const auto& cb : b) {
                if (ca.size() && cb.size() > std::numeric_limits<std::size_t>::max() / ca.size())
                    throw std::overflow_error("multichannel screening: count overflow");
                const std::size_t count = ca.size() * cb.size();
                if (count > (std::numeric_limits<std::size_t>::max() - bytes) / sizeof(GaussianPairIntermediate))
                    throw std::overflow_error("multichannel screening: size overflow");
                bytes += count * sizeof(GaussianPairIntermediate);
            }
        return bytes;
    };

    const auto contraction_start = Clock::now();
    for (std::size_t ii = 0; ii < np; ++ii) {
        for (const auto j : deactivate[ii]) reciprocal_active.erase(j);
        const std::size_t i = traversal[ii];
        const auto& x = pairs[i];
        LocalPreparedShellPair outer;
        bool ready = false;
        const auto evaluate = [&](std::size_t jj) {
            ++stats.candidate_visits;
            const std::size_t j = traversal[jj];
            const auto& y = pairs[j];
            const std::size_t outer_estimate = estimate(x);
            const std::size_t inner_estimate = i == j ? 0 : estimate(y);
            if (outer_estimate > memory_limit_bytes ||
                inner_estimate > memory_limit_bytes - outer_estimate)
                throw std::runtime_error("multichannel screening: Gaussian-pair memory limit exceeded");
            if (!ready) {
                outer = prepare_local_shell_pair(shells[x.a], shells[x.b]);
                ++stats.prepared_pair_instances;
                ready = true;
            }
            LocalPreparedShellPair inner;
            const LocalPreparedShellPair* second = &outer;
            if (i != j) {
                inner = prepare_local_shell_pair(shells[y.a], shells[y.b]);
                ++stats.prepared_pair_instances;
                second = &inner;
            }
            const std::size_t resident = pair_bytes(outer) + (i == j ? 0 : pair_bytes(inner));
            stats.peak_pair_bytes = std::max(stats.peak_pair_bytes, resident);
            if (resident > memory_limit_bytes)
                throw std::runtime_error("multichannel screening: Gaussian-pair memory limit exceeded");
            const auto quartet = local_shell_quartet_paired(outer, *second);
            ++stats.evaluated_quartets;
            for (std::size_t k : active)
                detail::accumulate_local_density_quartet(
                    out.matrices[k], densities[k], offsets,
                    x.a, x.b, y.a, y.b, i == j, quartet);
        };
        // Prefix is never directionally certified for destination ii.
        const std::size_t prefix_end = std::max(ii, cut[ii]);
        for (std::size_t jj = ii; jj < prefix_end; ++jj) evaluate(jj);
        // In the suffix, evaluate only ranks whose reciprocal test fails.
        for (auto it = reciprocal_active.lower_bound(prefix_end);
             it != reciprocal_active.end(); ++it) evaluate(*it);
    }
    stats.evaluated_quartets = stats.candidate_visits;
    stats.skipped_quartets = stats.total_quartets - stats.evaluated_quartets;
    stats.cutoff_skipped_quartets = stats.skipped_quartets;
    // For every output pair, mutual omissions are a subset of its certified
    // directed suffix. Its full suffix bound therefore certifies the error.
    for (std::size_t c = 0; c < active.size(); ++c) {
        const auto k = active[c];
        long double maximum = 0.0L;
        for (std::size_t r = 0; r < np; ++r)
            maximum = std::max(maximum, pairs[traversal[r]].bmax * suffix[c][cut[r]]);
        stats.max_error_bounds[k] = static_cast<double>(maximum);
    }
    for (std::size_t k : active)
        out.matrices[k] = Eigen::MatrixXd(out.matrices[k].selfadjointView<Eigen::Upper>());
    stats.contraction_seconds = std::chrono::duration<double>(Clock::now() - contraction_start).count();
    stats.total_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}

} // namespace chem::kumo
