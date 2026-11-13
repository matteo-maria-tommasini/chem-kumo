#pragma once

#include "chem/kumo/surface.hpp"

#include <Eigen/Core>
#include <functional>
#include <span>
#include <vector>

namespace chem::kumo::detail {

/// Values and gradients of a scalar function at the columns of `points`;
/// fills values[p] and gradients.col(p), already sized to the point count.
using ScalarFunction = std::function<void(const Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                                          std::vector<double>& values,
                                          Eigen::Matrix<double, 3, Eigen::Dynamic>& gradients)>;

/// Isosurface mesh with, for each vertex, the grid point (i, j, k) at the
/// lower end of the grid edge it lies on.
struct Isosurface {
    TriangleMesh mesh;
    std::vector<Eigen::Vector3i> vertex_cells;
};

/// Whether a sample on the outer faces of the grid is >= isovalue, so that
/// the region f >= isovalue may extend beyond the grid. The grid must be
/// validated and `samples` in Field order.
[[nodiscard]] bool isosurface_reaches_boundary(const Grid& grid, std::span<const double> samples,
                                               double isovalue);

/// Boundary of the region f >= isovalue, from samples of f on a grid (in
/// Field order) and the function itself.
///
/// Each grid cell is split into six tetrahedra around its main diagonal
/// (from corner (0,0,0) to (1,1,1)); these tile space conformingly, so
/// neighbouring cells share their faces and edges. Every tetrahedron edge
/// whose ends lie on either side of the isovalue (by the samples) carries
/// one vertex, shared by all tetrahedra around the edge; each tetrahedron
/// contributes a triangle (one corner apart) or a quadrilateral (two and
/// two), oriented outward (towards f < isovalue) and split along its
/// shorter diagonal. The vertex is the root of f - isovalue on its edge,
/// found by Newton steps on f and its gradient, safeguarded by bisection,
/// to 1e-10 bohr; the function is called on all unconverged vertices at
/// once. Vertices are numbered in order of first appearance, cells with k
/// fastest, so the mesh does not depend on threads.
///
/// Throws std::invalid_argument if the grid is invalid or has fewer than two
/// points along an axis, `samples` has the wrong size or a nonfinite value,
/// the isovalue is not finite, or isosurface_reaches_boundary() holds (the
/// surface would not be closed).
[[nodiscard]] Isosurface extract_isosurface(const Grid& grid, std::span<const double> samples,
                                            double isovalue, const ScalarFunction& function);

} // namespace chem::kumo::detail
