#include "boys_function.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace chem::kumo::detail {
namespace {
constexpr double table_spacing = 0.05;
constexpr double table_end = 40.0;
constexpr int taylor_terms = 7; // k = 0, ..., 6
// The Taylor series for F_m needs F_m, ..., F_(m+6) at the table point.
constexpr int table_orders = boys_max_order + taylor_terms;
constexpr int table_points = 801; // T = 0, 0.05, ..., 40

// F_m(T) for m = 0, ..., table_orders - 1 at every table point. The highest
// order comes from the series
//   F_M(T) = exp(-T) sum_k (2T)^k / ((2M + 1)(2M + 3) ... (2M + 2k + 1)),
// whose terms are all positive (no cancellation), and the lower orders from
// the stable downward recursion.
const std::vector<double>& table() {
    static const std::vector<double> values = [] {
        std::vector<double> result(static_cast<std::size_t>(table_points * table_orders));
        constexpr int top = table_orders - 1;
        for (int i = 0; i < table_points; ++i) {
            const double T = i * table_spacing;
            const double decay = std::exp(-T);
            double term = 1.0 / (2.0 * top + 1.0);
            double sum = term;
            for (int k = 1; k < 10000 && term > 1e-18 * sum; ++k) {
                term *= 2.0 * T / (2.0 * top + 2.0 * k + 1.0);
                sum += term;
            }
            double* row = result.data() + static_cast<std::size_t>(i * table_orders);
            row[top] = decay * sum;
            for (int m = top; m > 0; --m)
                row[m - 1] = (2.0 * T * row[m] + decay) / (2.0 * m - 1.0);
        }
        return result;
    }();
    return values;
}
} // namespace

void boys_function(double T, std::span<double> values) {
    if (values.empty() || values.size() > static_cast<std::size_t>(boys_max_order) + 1)
        throw std::invalid_argument("unsupported Boys function order");
    if (!std::isfinite(T) || T < 0.0)
        throw std::invalid_argument("Boys function argument must be finite and nonnegative");
    const int top = static_cast<int>(values.size()) - 1;
    const double decay = std::exp(-T);
    if (T >= table_end) {
        values[0] = 0.5 * std::sqrt(std::numbers::pi / T);
        for (int m = 0; m < top; ++m)
            values[static_cast<std::size_t>(m + 1)] =
                ((2.0 * m + 1.0) * values[static_cast<std::size_t>(m)] - decay) / (2.0 * T);
        return;
    }
    const int point = static_cast<int>(std::lround(T / table_spacing));
    const double delta = T - point * table_spacing; // |delta| <= 0.025
    const double* row = table().data() + static_cast<std::size_t>(point * table_orders);
    // Taylor expansion about the nearest table point T0:
    // F_m(T0 + delta) = sum_k F_(m+k)(T0) (-delta)^k / k!, since dF_m/dT = -F_(m+1).
    double sum = 0.0;
    double factor = 1.0;
    for (int k = 0; k < taylor_terms; ++k) {
        sum += row[top + k] * factor;
        factor *= -delta / (k + 1.0);
    }
    values[static_cast<std::size_t>(top)] = sum;
    for (int m = top; m > 0; --m)
        values[static_cast<std::size_t>(m - 1)] =
            (2.0 * T * values[static_cast<std::size_t>(m)] + decay) / (2.0 * m - 1.0);
}
} // namespace chem::kumo::detail
