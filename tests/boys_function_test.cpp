// Boys function F_m(T) = integral_0^1 t^(2m) exp(-T t^2) dt against
// independent references: Romberg quadrature of the definition, the closed
// forms F_m(0) = 1 / (2m + 1) and F_0(T) = sqrt(pi / T) erf(sqrt(T)) / 2, and
// the recurrence F_(m+1) = ((2m + 1) F_m - exp(-T)) / (2T).
#include "boys_function.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using chem::kumo::detail::boys_function;
using chem::kumo::detail::boys_max_order;

// Romberg integration of t^(2m) exp(-T t^2) over [0, 1]: trapezoid rules with
// 1, 2, 4, ..., 2^20 intervals and Richardson extrapolation, stopped when two
// successive diagonal entries agree to 1e-16 relative.
double romberg(int m, double T) {
    const auto f = [&](double t) { return std::pow(t, 2 * m) * std::exp(-T * t * t); };
    std::vector<double> previous, current;
    double trapezoid = 0.5 * (f(0.0) + f(1.0));
    previous.push_back(trapezoid);
    for (int level = 1; level <= 20; ++level) {
        const long intervals = 1L << level;
        const double h = 1.0 / static_cast<double>(intervals);
        // Compensated (Kahan) sum: up to 2^19 terms would otherwise lose
        // about 1e-13 to rounding.
        double midpoints = 0.0, compensation = 0.0;
        for (long i = 1; i < intervals; i += 2) {
            const double y = f(static_cast<double>(i) * h) - compensation;
            const double t = midpoints + y;
            compensation = (t - midpoints) - y;
            midpoints = t;
        }
        trapezoid = 0.5 * trapezoid + h * midpoints;
        current.assign(1, trapezoid);
        double scale = 1.0;
        for (int k = 1; k <= level; ++k) {
            scale *= 4.0;
            current.push_back(current[static_cast<std::size_t>(k - 1)] +
                              (current[static_cast<std::size_t>(k - 1)] -
                               previous[static_cast<std::size_t>(k - 1)]) /
                                  (scale - 1.0));
        }
        if (level >= 6 &&
            std::abs(current.back() - previous.back()) <= 1e-16 * std::abs(current.back()))
            return current.back();
        previous = current;
    }
    return previous.back();
}

double relative(double actual, double expected) {
    return std::abs(actual - expected) / std::max(std::abs(expected), 1e-300);
}
} // namespace

int main() {
    std::array<double, boys_max_order + 1> values{};

    // F_m(0) = 1 / (2m + 1).
    boys_function(0.0, values);
    for (int m = 0; m <= boys_max_order; ++m)
        kumo_test::check(relative(values[static_cast<std::size_t>(m)], 1.0 / (2.0 * m + 1.0)) <
                             1e-15,
                         "F_" + std::to_string(m) + "(0) = 1/(2m+1)");

    // Against quadrature, including both sides of every switch: the Taylor
    // midpoints (0.025 from a table point) and the change of method at T = 40.
    const double arguments[] = {1e-14,  1e-8,  1e-3,   0.024, 0.025, 0.026, 0.3,  1.0,
                                2.5,    4.975, 7.3,    12.0,  19.99, 25.0,  33.3, 39.975,
                                39.999, 40.0,  40.001, 47.5,  60.0,  90.0,  120.0};
    double worst = 0.0;
    for (const double T : arguments) {
        boys_function(T, values);
        for (int m = 0; m <= boys_max_order; ++m) {
            const double error = relative(values[static_cast<std::size_t>(m)], romberg(m, T));
            worst = std::max(worst, error);
            kumo_test::check(error < 1e-14, "F_" + std::to_string(m) + "(" + kumo_test::num(T) +
                                                ") agrees with quadrature (relative error " +
                                                kumo_test::num(error) + ")");
        }
    }
    std::cout << "largest relative difference from quadrature: " << worst << '\n';

    // F_0 in closed form over a fine scan, and the upward recurrence where it
    // is well conditioned (2T > 2m + 1).
    for (double T = 1e-6; T < 400.0; T = T * 1.07 + 0.003) {
        boys_function(T, values);
        const double f0 = 0.5 * std::sqrt(std::numbers::pi / T) * std::erf(std::sqrt(T));
        kumo_test::check(relative(values[0], f0) < 1e-14,
                         "F_0(" + kumo_test::num(T) + ") = sqrt(pi/T) erf(sqrt(T)) / 2");
        for (int m = 0; m < boys_max_order && 2.0 * T > 4.0 * m + 2.0; ++m) {
            const double upward =
                ((2.0 * m + 1.0) * values[static_cast<std::size_t>(m)] - std::exp(-T)) / (2.0 * T);
            kumo_test::check(relative(values[static_cast<std::size_t>(m + 1)], upward) < 1e-13,
                             "F_" + std::to_string(m + 1) + "(" + kumo_test::num(T) +
                                 ") satisfies the recurrence");
        }
    }

    // Fewer orders give the same leading values.
    boys_function(3.7, values);
    std::array<double, 4> few{};
    boys_function(3.7, few);
    for (std::size_t m = 0; m < few.size(); ++m)
        kumo_test::check(relative(few[m], values[m]) < 1e-14,
                         "a shorter request gives the same F_" + std::to_string(m));

    // Invalid input fails.
    std::array<double, boys_max_order + 2> too_many{};
    kumo_test::check_throws<std::invalid_argument>([&] { boys_function(1.0, too_many); },
                                                   "too many orders are rejected");
    kumo_test::check_throws<std::invalid_argument>([&] { boys_function(1.0, std::span<double>{}); },
                                                   "an empty request is rejected");
    kumo_test::check_throws<std::invalid_argument>([&] { boys_function(-1e-12, values); },
                                                   "a negative argument is rejected");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { boys_function(std::numeric_limits<double>::quiet_NaN(), values); },
        "NaN is rejected");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { boys_function(std::numeric_limits<double>::infinity(), values); },
        "infinity is rejected");
    return kumo_test::exit_code();
}
