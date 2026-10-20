#pragma once

#include "chem/kumo/local_density_overlap_multichannel.hpp"
#include "chem/kumo/local_pair_bounds.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

struct LocalMultichannelScreeningStatistics {
    std::size_t shell_pairs = 0;
    std::size_t total_quartets = 0;
    std::size_t evaluated_quartets = 0;
    std::size_t skipped_quartets = 0;
    std::size_t prepared_pair_instances = 0;
    std::size_t peak_pair_bytes = 0;
    std::size_t active_channels = 0;
    std::vector<double> max_error_bounds;
    double norm_seconds = 0.0;
    double contraction_seconds = 0.0;
    double total_seconds = 0.0;
};

struct LocalMultichannelScreenedDensityOverlap {
    std::vector<Eigen::MatrixXd> matrices;
    LocalMultichannelScreeningStatistics statistics;
};

// Canonical streaming traversal. A quartet is skipped only when every active
// channel can charge both affected shell-pair error budgets. The memory cap
// applies to prepared Gaussian-pair vectors, not to total process RSS.
[[nodiscard]] inline LocalMultichannelScreenedDensityOverlap
local_density_overlap_multichannel_screened(
    const BasisSet& basis,
    const std::vector<Eigen::MatrixXd>& densities,
    const std::vector<double>& epsilons,
    std::size_t memory_limit_bytes) {
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
    LocalMultichannelScreenedDensityOverlap out;
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
    std::vector<std::vector<long double>> spent(densities.size(), std::vector<long double>(np, 0.0L));

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
    for (std::size_t i = 0; i < np; ++i) {
        const auto& x = pairs[i];
        const long double mi = x.a == x.b ? 1.0L : 2.0L;
        LocalPreparedShellPair outer;
        bool ready = false;
        for (std::size_t j = i; j < np; ++j) {
            const auto& y = pairs[j];
            const long double mj = y.a == y.b ? 1.0L : 2.0L;
            bool skip = true;
            for (std::size_t k : active) {
                const long double to_i = mj * x.bmax * dsum[k][j];
                const long double to_j = i == j ? 0.0L : mi * y.bmax * dsum[k][i];
                if (spent[k][i] + to_i > static_cast<long double>(epsilons[k]) ||
                    (i != j && spent[k][j] + to_j > static_cast<long double>(epsilons[k]))) {
                    skip = false;
                    break;
                }
            }
            if (skip) {
                for (std::size_t k : active) {
                    spent[k][i] += mj * x.bmax * dsum[k][j];
                    if (i != j) spent[k][j] += mi * y.bmax * dsum[k][i];
                }
                ++stats.skipped_quartets;
                continue;
            }
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
        }
    }
    for (std::size_t k : active)
        out.matrices[k] = Eigen::MatrixXd(out.matrices[k].selfadjointView<Eigen::Upper>());
    for (std::size_t k : active)
        stats.max_error_bounds[k] = static_cast<double>(*std::max_element(spent[k].begin(), spent[k].end()));
    stats.contraction_seconds = std::chrono::duration<double>(Clock::now() - contraction_start).count();
    stats.total_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}

} // namespace chem::kumo
