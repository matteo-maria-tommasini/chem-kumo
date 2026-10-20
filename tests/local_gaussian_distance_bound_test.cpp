#include "chem/kumo/local_gaussian_distance_bound.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

int main() {
    using namespace chem::kumo;
    std::mt19937_64 rng(0xC010CULL);
    std::uniform_real_distribution<double> pos(-5., 5.), expo(-2., 2.), scale(-1., 1.);
    std::uniform_int_distribution<int> degree(0, 3);
    double max_violation = 0., min_ratio = 1e300, max_ratio = 0.;
    int checked = 0, separated = 0;
    for (int k = 0; k < 5000; ++k) {
        std::array<LocalGaussianPrimitive, 4> g{};
        for (auto& x : g) {
            x.exponent = std::exp(expo(rng));
            x.scale = scale(rng);
            for (int d = 0; d < 3; ++d) {
                x.center[d] = pos(rng);
                x.powers[d] = degree(rng);
            }
        }
        const auto a = make_gaussian_pair(g[0], g[1]);
        const auto b = make_gaussian_pair(g[2], g[3]);
        const double exact = local_gaussian_pair_integral(a, b);
        const double upper = local_gaussian_pair_distance_bound(a, b);
        if (!std::isfinite(exact) || !std::isfinite(upper) || upper < 0)
            throw std::runtime_error("nonfinite integral/bound");
        const double violation = std::max(0., std::abs(exact) - upper);
        max_violation = std::max(max_violation, violation);
        if (violation > 5e-12 * std::max({1., std::abs(exact), upper}))
            throw std::runtime_error("distance bound violated");
        if (upper > 1e-200 && std::abs(exact) > 1e-200) {
            const double ratio = upper / std::abs(exact);
            min_ratio = std::min(min_ratio, ratio);
            max_ratio = std::max(max_ratio, ratio);
            ++checked;
        }
        if (std::hypot(a.center[0]-b.center[0], a.center[1]-b.center[1]) > 3)
            ++separated;
    }
    std::cout << std::setprecision(12) << "cases=5000 checked_ratios=" << checked
              << " separated=" << separated << " max_positive_violation=" << max_violation
              << " min_bound_to_integral=" << min_ratio
              << " max_bound_to_integral=" << max_ratio << '\n';
}
