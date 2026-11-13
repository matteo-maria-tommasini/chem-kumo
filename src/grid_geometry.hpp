#pragma once

#include "chem/kumo/density.hpp"

#include <Eigen/LU> // Matrix3d::determinant()
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace chem::kumo::detail {

/// Validate a regular grid and return its number of points.
///
/// This is the single geometry rule for every library entry point that
/// evaluates on a Grid. It requires:
///   - positive point counts whose product fits in std::size_t;
///   - a finite origin and finite step vectors;
///   - nonzero step vectors;
///   - step directions that are not (nearly) linearly dependent.
///
/// The independence test uses the determinant of the *normalized* step
/// directions, so it depends only on the angles between the axes and not on
/// the spacing: a grid with 1e-7 bohr spacing is as valid as one with 0.2 bohr.
/// Throws std::invalid_argument on any violation.
///
/// Format-specific limits (for example the minimum of two points per axis
/// required by Cube files or trilinear sampling) are checked by their callers.
inline std::size_t validated_point_count(const Grid& grid) {
    std::size_t count = 1;
    for (int axis = 0; axis < 3; ++axis) {
        const int n = grid.counts[axis];
        if (n <= 0 || static_cast<std::size_t>(n) > std::numeric_limits<std::size_t>::max() / count)
            throw std::invalid_argument("invalid or oversized grid");
        count *= static_cast<std::size_t>(n);
    }
    if (!grid.origin.allFinite() || !grid.steps.allFinite())
        throw std::invalid_argument("invalid grid geometry: nonfinite origin or step");

    Eigen::Matrix3d directions;
    for (int axis = 0; axis < 3; ++axis) {
        const double length = grid.steps.col(axis).norm();
        if (!std::isfinite(length) || length <= 0.0)
            throw std::invalid_argument("invalid grid geometry: zero step vector");
        directions.col(axis) = grid.steps.col(axis) / length;
    }
    // |det| of unit vectors is 1 for orthogonal axes and 0 for dependent ones.
    if (std::abs(directions.determinant()) < 1e-12)
        throw std::invalid_argument("invalid grid geometry: linearly dependent axes");
    return count;
}

/// Cartesian position (bohr) of grid point number `linear`.
///
/// Points are numbered with k (the third axis) fastest, then j, then i:
/// linear = (i * ny + j) * nz + k. This is the storage order of Field values
/// and of Cube files. The grid must already be validated.
inline Eigen::Vector3d grid_point(const Grid& grid, std::size_t linear) {
    const auto ny = static_cast<std::size_t>(grid.counts[1]);
    const auto nz = static_cast<std::size_t>(grid.counts[2]);
    const std::size_t i = linear / (ny * nz);
    const std::size_t j = (linear / nz) % ny;
    const std::size_t k = linear % nz;
    return grid.origin + grid.steps * Eigen::Vector3d(static_cast<double>(i),
                                                      static_cast<double>(j),
                                                      static_cast<double>(k));
}

/// Linear index of grid point (i, j, k) in Field order (k fastest).
inline std::size_t grid_index(const Grid& grid, int i, int j, int k) {
    return (static_cast<std::size_t>(i) * static_cast<std::size_t>(grid.counts[1]) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(grid.counts[2]) +
           static_cast<std::size_t>(k);
}

/// Cartesian positions (bohr) of grid points start, ..., start + count - 1,
/// one per column.
inline Eigen::Matrix<double, 3, Eigen::Dynamic> grid_points(const Grid& grid, std::size_t start,
                                                            std::size_t count) {
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(count));
    for (std::size_t local = 0; local < count; ++local)
        points.col(static_cast<Eigen::Index>(local)) = grid_point(grid, start + local);
    return points;
}

} // namespace chem::kumo::detail
