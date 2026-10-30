#pragma once

#include "chem/kumo/density.hpp"

#include <cmath>
#include <stdexcept>

namespace chem::kumo {

/// Volume of a cell in the affine grid (bohr^3).
[[nodiscard]] inline double grid_cell_volume(const Field& field) {
    const double volume = std::abs(field.grid().steps.determinant());
    if (!(volume > 0.0) || !std::isfinite(volume))
        throw std::invalid_argument("unrepresentable grid cell volume");
    return volume;
}

/// Integrate a scalar field with a rectangular rule, including all stored samples.
[[nodiscard]] inline double integrate_rectangular(const Field& field) {
    const double volume = grid_cell_volume(field);
    long double sum = 0.0L;
    for (const double value : field.values()) {
        if (!std::isfinite(value))
            throw std::invalid_argument("nonfinite field sample");
        sum += value;
    }
    const double result = static_cast<double>(sum * volume);
    if (!std::isfinite(result))
        throw std::overflow_error("field integral overflow");
    return result;
}

/// Integrate a scalar field with tensor-product trapezoidal weights.
/// A singleton grid axis has unit weight (no interpolation is required).
[[nodiscard]] inline double integrate_trapezoidal(const Field& field) {
    const double volume = grid_cell_volume(field);
    const Eigen::Vector3i counts = field.grid().counts;
    long double sum = 0.0L;
    for (int i = 0; i < counts[0]; ++i)
        for (int j = 0; j < counts[1]; ++j)
            for (int k = 0; k < counts[2]; ++k) {
                const double value = field(i, j, k);
                if (!std::isfinite(value))
                    throw std::invalid_argument("nonfinite field sample");
                const double wi = counts[0] == 1 ? 1.0 : ((i == 0 || i == counts[0]-1) ? 0.5 : 1.0);
                const double wj = counts[1] == 1 ? 1.0 : ((j == 0 || j == counts[1]-1) ? 0.5 : 1.0);
                const double wk = counts[2] == 1 ? 1.0 : ((k == 0 || k == counts[2]-1) ? 0.5 : 1.0);
                sum += wi * wj * wk * value;
            }
    const double result = static_cast<double>(sum * volume);
    if (!std::isfinite(result))
        throw std::overflow_error("field integral overflow");
    return result;
}

} // namespace chem::kumo
