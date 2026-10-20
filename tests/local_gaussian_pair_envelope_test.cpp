#include "chem/kumo/local_gaussian_pair_envelope.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

int main() {
    using namespace chem::kumo;
    std::mt19937_64 rng(0xC010EULL);
    std::uniform_real_distribution<double> position(-4., 4.), exponent(-1.5, 1.5), scale(-1., 1.);
    std::uniform_int_distribution<int> angular(0, 3);
    constexpr std::array<double, 4> gammas{0.25, 0.5, 0.7, 0.85};
    std::array<int, 4> tighter{};
    double worst_point_violation = 0., worst_integral_violation = 0.;
    for (int k = 0; k < 3000; ++k) {
        std::array<LocalGaussianPrimitive, 4> primitives{};
        for (auto& g : primitives) {
            g.exponent = std::exp(exponent(rng));
            g.scale = scale(rng);
            for (int d = 0; d < 3; ++d) {
                g.center[d] = position(rng);
                g.powers[d] = angular(rng);
            }
        }
        const auto a = make_gaussian_pair(primitives[0], primitives[1]);
        const auto b = make_gaussian_pair(primitives[2], primitives[3]);
        const double integral = std::abs(local_gaussian_pair_integral(a, b));
        for (std::size_t j = 0; j < gammas.size(); ++j) {
            const auto ea = make_gaussian_pair_envelope(a, gammas[j]);
            const auto eb = make_gaussian_pair_envelope(b, gammas[j]);
            const double bound = local_gaussian_envelope_overlap_bound(ea, eb);
            worst_integral_violation = std::max(worst_integral_violation, integral - bound);
            if (bound < 1e-12) ++tighter[j];
            for (const auto* pair : {&a, &b}) {
                const auto& env = pair == &a ? ea : eb;
                std::array<double, 3> r{};
                for (int d = 0; d < 3; ++d) r[d] = position(rng);
                double polynomial = pair->prefactor, r2 = 0., env_r2 = 0.;
                for (int d = 0; d < 3; ++d) {
                    const double x = r[d] - pair->center[d];
                    double axis = 0.;
                    for (int n = pair->degree[d]; n >= 0; --n)
                        axis = axis * x + pair->polynomial[d][n];
                    polynomial *= axis;
                    r2 += x * x;
                    const double y = r[d] - env.center[d];
                    env_r2 += y * y;
                }
                const double actual = std::abs(polynomial * std::exp(-pair->exponent * r2));
                const double upper = env.amplitude * std::exp(-env.exponent * env_r2);
                worst_point_violation = std::max(worst_point_violation, actual - upper);
            }
        }
    }
    std::cout << std::setprecision(12) << "cases=3000 max_point_violation=" << worst_point_violation
              << " max_integral_violation=" << worst_integral_violation;
    for (std::size_t j = 0; j < gammas.size(); ++j)
        std::cout << " gamma=" << gammas[j] << " below_1e-12=" << tighter[j];
    std::cout << '\n';
    if (worst_point_violation > 1e-10 || worst_integral_violation > 1e-10)
        throw std::runtime_error("Gaussian envelope violation");
}
