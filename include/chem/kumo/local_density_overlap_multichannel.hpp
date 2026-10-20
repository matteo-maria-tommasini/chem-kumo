#pragma once

#include "chem/kumo/local_density_quartet_contraction.hpp"
#include "chem/kumo/local_shell_pair_integrals.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

struct LocalMultichannelStatistics {
    std::size_t shell_pairs = 0;
    std::size_t active_channels = 0;
    std::size_t evaluated_quartets = 0;
    std::size_t prepared_pair_instances = 0;
    double quartet_seconds = 0.0;
    double accumulation_seconds = 0.0;
    double total_seconds = 0.0;
};

struct LocalMultichannelDensityOverlap {
    std::vector<Eigen::MatrixXd> matrices;
    LocalMultichannelStatistics statistics;
};

// Unscreened reference implementation: one quartet evaluation shared by all channels.
// Pair storage is bounded to the outer and current inner shell pair, independently
// of the number of shell quartets. Densities must be real, finite and symmetric.
[[nodiscard]] inline LocalMultichannelDensityOverlap local_density_overlap_multichannel(
    const BasisSet& basis, const std::vector<Eigen::MatrixXd>& densities) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao = offsets.back();
    for (const auto& p : densities) {
        if (p.rows() != nao || p.cols() != nao || !p.allFinite())
            throw std::invalid_argument("multichannel overlap: invalid density");
        if (nao && (p - p.transpose()).cwiseAbs().maxCoeff() >
                        1e-10 * std::max(1.0, p.cwiseAbs().maxCoeff()))
            throw std::invalid_argument("multichannel overlap: asymmetric density");
    }
    LocalMultichannelDensityOverlap out;
    out.matrices.reserve(densities.size());
    for (std::size_t k = 0; k < densities.size(); ++k)
        out.matrices.emplace_back(Eigen::MatrixXd::Zero(nao, nao));
    // Exact-zero detection only: never discard a small but nonzero density.
    // Validation above must run even for channels that turn out to be inactive.
    std::vector<std::size_t> active_channels;
    active_channels.reserve(densities.size());
    for (std::size_t k = 0; k < densities.size(); ++k)
        if (!(densities[k].array() == 0.0).all())
            active_channels.push_back(k);
    out.statistics.active_channels = active_channels.size();
    if (active_channels.empty() || nao == 0) {
        out.statistics.total_seconds = std::chrono::duration<double>(Clock::now() - start).count();
        return out;
    }
    struct Pair { std::size_t a, b; };
    std::vector<Pair> pairs;
    for (std::size_t a = 0; a < shells.size(); ++a)
        for (std::size_t b = a; b < shells.size(); ++b)
            pairs.push_back({a, b});
    out.statistics.shell_pairs = pairs.size();
    // Keep the established canonical traversal and quartet accumulation order.
    for (std::size_t ip = 0; ip < pairs.size(); ++ip) {
        const auto [a, b] = pairs[ip];
        const auto outer = prepare_local_shell_pair(shells[a], shells[b]);
        ++out.statistics.prepared_pair_instances;
        for (std::size_t iq = ip; iq < pairs.size(); ++iq) {
            const auto [c, d] = pairs[iq];
            LocalPreparedShellPair inner;
            const LocalPreparedShellPair* second = &outer;
            if (ip != iq) {
                inner = prepare_local_shell_pair(shells[c], shells[d]);
                ++out.statistics.prepared_pair_instances;
                second = &inner;
            }
            const auto tq = Clock::now();
            const auto quartet = local_shell_quartet_paired(outer, *second);
            out.statistics.quartet_seconds += std::chrono::duration<double>(Clock::now() - tq).count();
            ++out.statistics.evaluated_quartets;
            const auto ta = Clock::now();
            for (const std::size_t channel : active_channels)
                detail::accumulate_local_density_quartet(
                    out.matrices[channel], densities[channel], offsets,
                    a, b, c, d, ip == iq, quartet);
            out.statistics.accumulation_seconds += std::chrono::duration<double>(Clock::now() - ta).count();
        }
    }
    for (auto& matrix : out.matrices)
        matrix = Eigen::MatrixXd(matrix.selfadjointView<Eigen::Upper>());
    out.statistics.total_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}

} // namespace chem::kumo
