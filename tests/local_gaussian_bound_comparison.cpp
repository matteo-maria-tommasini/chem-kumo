#include "chem/kumo/local_gaussian_distance_bound.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
using Clock = std::chrono::steady_clock;

struct Case {
    GaussianPairIntermediate a, b;
    double cs, exact;
    bool separated;
};

int main() {
    constexpr int n = 6000;
    std::mt19937_64 rng(0xC010DULL);
    std::uniform_real_distribution<double> near(-1., 1.), far(-7., 7.);
    std::uniform_real_distribution<double> exponent(-1.5, 1.5), scale(0.2, 1.);
    std::uniform_int_distribution<int> power(0, 2);
    std::vector<Case> cases;
    cases.reserve(n);
    for (int k = 0; k < n; ++k) {
        std::array<LocalGaussianPrimitive, 4> g{};
        const bool separated = (k % 2 == 1);
        std::array<double, 3> origin{};
        for (int d = 0; d < 3; ++d) origin[d] = near(rng);
        for (int j = 0; j < 4; ++j) {
            auto& v = g[j];
            v.exponent = std::exp(exponent(rng));
            v.scale = scale(rng) * ((k + j) % 3 == 0 ? -1. : 1.);
            for (int d = 0; d < 3; ++d) {
                v.center[d] = origin[d] + 0.25 * near(rng)
                    + (separated && j >= 2 ? far(rng) : 0.);
                v.powers[d] = power(rng);
            }
        }
        auto a = make_gaussian_pair(g[0], g[1]);
        auto b = make_gaussian_pair(g[2], g[3]);
        const double aa = local_gaussian_pair_integral(a, a);
        const double bb = local_gaussian_pair_integral(b, b);
        if (aa < -1e-12 || bb < -1e-12) throw std::runtime_error("negative pair norm");
        const double cs = std::sqrt(std::max(0., aa)) * std::sqrt(std::max(0., bb));
        const double exact = local_gaussian_pair_integral(a, b);
        cases.push_back({a, b, cs, exact, separated});
    }
    std::size_t tighter = 0, near_tighter = 0, far_tighter = 0;
    std::size_t near_count = 0, far_count = 0;
    double max_violation = 0., worst_ratio = 0.;
    double sum_log_ratio = 0.;
    std::size_t valid_ratios = 0;
    std::array<std::size_t, 4> cs_skips{}, combined_skips{};
    const std::array<double, 4> thresholds{1e-6, 1e-8, 1e-10, 1e-12};
    for (const auto& c : cases) {
        const double db = local_gaussian_pair_distance_bound(c.a, c.b);
        if (!std::isfinite(db) || db < 0) throw std::runtime_error("invalid distance bound");
        const double violation = std::max(0., std::abs(c.exact) - db);
        max_violation = std::max(max_violation, violation);
        if (violation > 1e-11 * std::max({1., db, std::abs(c.exact)}))
            throw std::runtime_error("distance bound violated");
        if (std::abs(c.exact) > c.cs + 1e-11 * std::max(1., c.cs))
            throw std::runtime_error("CS bound violated");
        const bool better = db < c.cs;
        tighter += better;
        if (c.separated) { ++far_count; far_tighter += better; }
        else { ++near_count; near_tighter += better; }
        if (c.cs > 1e-250 && db > 0) {
            const double ratio = db / c.cs;
            worst_ratio = std::max(worst_ratio, ratio);
            sum_log_ratio += std::log(ratio);
            ++valid_ratios;
        }
        for (std::size_t j = 0; j < thresholds.size(); ++j) {
            cs_skips[j] += c.cs < thresholds[j];
            combined_skips[j] += std::min(c.cs, db) < thresholds[j];
        }
    }
    volatile double sink = 0.;
    auto time_kernel = [&](bool bound) {
        const auto start = Clock::now();
        double checksum = 0.;
        for (int repeat = 0; repeat < 3; ++repeat)
            for (const auto& c : cases)
                checksum += bound ? local_gaussian_pair_distance_bound(c.a, c.b)
                                  : local_gaussian_pair_integral(c.a, c.b);
        sink = checksum;
        return std::chrono::duration<double>(Clock::now() - start).count() / 3.;
    };
    const double integral_seconds = time_kernel(false);
    const double bound_seconds = time_kernel(true);
    std::cout << std::setprecision(12)
              << "cases=" << cases.size() << " near=" << near_count
              << " separated=" << far_count << " tighter=" << tighter
              << " near_tighter=" << near_tighter << " separated_tighter=" << far_tighter
              << " max_positive_violation=" << max_violation
              << " geometric_mean_distance_to_cs=" << std::exp(sum_log_ratio / valid_ratios)
              << " max_distance_to_cs=" << worst_ratio << '\n';
    for (std::size_t j = 0; j < thresholds.size(); ++j)
        std::cout << "threshold=" << thresholds[j]
                  << " cs_skips=" << cs_skips[j]
                  << " combined_skips=" << combined_skips[j] << '\n';
    std::cout << "integral_seconds=" << integral_seconds
              << " distance_bound_seconds=" << bound_seconds
              << " cost_ratio=" << bound_seconds / integral_seconds
              << " checksum=" << sink << '\n';
}
