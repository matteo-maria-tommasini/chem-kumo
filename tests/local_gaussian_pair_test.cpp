#include "chem/kumo/local_gaussian_pair.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

int main() {
    using namespace chem::kumo;
    std::mt19937_64 rng(0xC0939BULL);
    std::uniform_real_distribution<double> xyz(-2.0, 2.0);
    std::uniform_real_distribution<double> expo(-2.0, 2.0);
    std::uniform_real_distribution<double> scale(-1.0, 1.0);
    std::uniform_int_distribution<int> angular(0, 3);
    constexpr int count = 4000;
    std::vector<std::array<LocalGaussianPrimitive, 4>> samples;
    samples.reserve(count);
    for (int i = 0; i < count; ++i) {
        std::array<LocalGaussianPrimitive, 4> quartet;
        for (auto& g : quartet) {
            g.exponent = std::exp(expo(rng));
            g.scale = scale(rng);
            for (int d = 0; d < 3; ++d) {
                g.center[d] = xyz(rng);
                g.powers[d] = angular(rng);
            }
        }
        samples.push_back(quartet);
    }
    std::vector<std::array<GaussianPairIntermediate, 2>> pairs;
    pairs.reserve(samples.size());
    double max_abs = 0.0, max_scaled = 0.0;
    for (const auto& q : samples) {
        const auto a = make_gaussian_pair(q[0], q[1]);
        const auto b = make_gaussian_pair(q[2], q[3]);
        const double reference = local_four_gaussian_integral(q);
        const double actual = local_gaussian_pair_integral(a, b);
        const double diff = std::abs(reference - actual);
        max_abs = std::max(max_abs, diff);
        max_scaled = std::max(max_scaled, diff / std::max({1.0, std::abs(reference), std::abs(actual)}));
        if (diff > 1e-11 * std::max({1.0, std::abs(reference), std::abs(actual)}))
            throw std::runtime_error("pair kernel disagrees with direct kernel");
        pairs.push_back({a, b});
    }
    volatile double checksum = 0.0;
    constexpr int repeats = 5;
    const auto begin_direct = std::chrono::steady_clock::now();
    for (int repeat = 0; repeat < repeats; ++repeat)
        for (const auto& q : samples) checksum = checksum + local_four_gaussian_integral(q);
    const auto end_direct = std::chrono::steady_clock::now();
    const auto begin_pair = std::chrono::steady_clock::now();
    for (int repeat = 0; repeat < repeats; ++repeat)
        for (const auto& p : pairs) checksum = checksum + local_gaussian_pair_integral(p[0], p[1]);
    const auto end_pair = std::chrono::steady_clock::now();
    const double direct_seconds = std::chrono::duration<double>(end_direct - begin_direct).count();
    const double pair_seconds = std::chrono::duration<double>(end_pair - begin_pair).count();
    std::cout << std::setprecision(12) << "cases=" << count
              << " max_abs=" << max_abs << " max_scaled=" << max_scaled
              << " direct_seconds=" << direct_seconds
              << " pair_seconds=" << pair_seconds
              << " speedup=" << direct_seconds / pair_seconds
              << " pair_bytes=" << sizeof(GaussianPairIntermediate)
              << " checksum=" << checksum << '\n';
}
