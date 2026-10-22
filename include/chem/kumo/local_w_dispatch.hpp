#pragma once

#include "chem/kumo/local_density_overlap_multichannel_hierarchical.hpp"

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chem::kumo {

enum class LocalWBackend { Reference, Screened, Hierarchical };

struct LocalWOptions {
    LocalWBackend backend = LocalWBackend::Reference;
    std::vector<double> error_budgets;
    std::size_t pair_memory_limit_bytes = 1024 * 1024;
    MultichannelScreeningOrder screening_order = MultichannelScreeningOrder::Canonical;
};

struct LocalWStatistics {
    std::size_t shell_pairs = 0;
    std::size_t total_quartets = 0;
    std::size_t evaluated_quartets = 0;
    std::size_t skipped_quartets = 0;
    std::size_t active_channels = 0;
    std::vector<double> max_error_bounds;
    double total_seconds = 0.0;
};

struct LocalWResult {
    std::vector<Eigen::MatrixXd> matrices;
    LocalWBackend backend = LocalWBackend::Reference;
    LocalWStatistics statistics;
};

// Thin, opt-in adapter: numerical algorithms and their legacy APIs are unchanged.
// The error bounds are those reported by the selected existing backend; they
// are not newly certified by this adapter.
[[nodiscard]] inline LocalWResult local_w_evaluate(
    const BasisSet& basis,
    const std::vector<Eigen::MatrixXd>& densities,
    const LocalWOptions& options = {}) {
    if (options.backend == LocalWBackend::Reference) {
        if (!options.error_budgets.empty())
            throw std::invalid_argument("local W: reference backend does not accept error budgets");
        auto result = local_density_overlap_multichannel(basis, densities);
        const auto& s = result.statistics;
        LocalWResult out;
        out.matrices = std::move(result.matrices);
        out.backend = options.backend;
        out.statistics.shell_pairs = s.shell_pairs;
        out.statistics.total_quartets = s.evaluated_quartets;
        out.statistics.evaluated_quartets = s.evaluated_quartets;
        out.statistics.active_channels = s.active_channels;
        out.statistics.max_error_bounds.assign(densities.size(), 0.0);
        out.statistics.total_seconds = s.total_seconds;
        return out;
    }
    if (options.backend != LocalWBackend::Screened &&
        options.backend != LocalWBackend::Hierarchical)
        throw std::invalid_argument("local W: unknown backend");
    if (options.error_budgets.size() != densities.size())
        throw std::invalid_argument("local W: error budget count mismatch");
    for (double e : options.error_budgets)
        if (!(e >= 0.0) || !std::isfinite(e))
            throw std::invalid_argument("local W: invalid error budget");
    if (options.backend == LocalWBackend::Screened) {
        auto result = local_density_overlap_multichannel_screened(
            basis, densities, options.error_budgets,
            options.pair_memory_limit_bytes, options.screening_order);
        const auto& s = result.statistics;
        LocalWResult out;
        out.matrices = std::move(result.matrices);
        out.backend = options.backend;
        out.statistics = {s.shell_pairs, s.total_quartets, s.evaluated_quartets,
                          s.skipped_quartets, s.active_channels,
                          s.max_error_bounds, s.total_seconds};
        return out;
    }
    auto result = local_density_overlap_multichannel_hierarchical(
        basis, densities, options.error_budgets,
        options.pair_memory_limit_bytes, options.screening_order);
    const auto& s = result.statistics;
    LocalWResult out;
    out.matrices = std::move(result.matrices);
    out.backend = options.backend;
    out.statistics = {s.shell_pairs, s.total_quartets, s.evaluated_quartets,
                      s.skipped_quartets, s.active_channels,
                      s.max_error_bounds, s.total_seconds};
    return out;
}

} // namespace chem::kumo
