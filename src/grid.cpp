// Grid constructors (around atoms, along a line, on a plane) and Field.

#include "chem/kumo/density.hpp"
#include "grid_geometry.hpp"

#include <Eigen/Geometry> // cross()
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace chem::kumo {
Grid Grid::around(const std::vector<Atom>& atoms, double margin, double spacing) {
    if (atoms.empty() || !std::isfinite(margin) || margin < 0 || !std::isfinite(spacing) ||
        spacing <= 0)
        throw std::invalid_argument("invalid grid parameters");
    if (!atoms.front().position.allFinite() || atoms.front().atomic_number <= 0)
        throw std::invalid_argument("invalid atom in grid");
    Eigen::Vector3d lo = atoms.front().position, hi = lo;
    for (const auto& a : atoms) {
        if (!a.position.allFinite() || a.atomic_number <= 0)
            throw std::invalid_argument("invalid atom in grid");
        lo = lo.cwiseMin(a.position);
        hi = hi.cwiseMax(a.position);
    }
    Grid g;
    g.origin = lo.array() - margin;
    for (int k = 0; k < 3; ++k) {
        const double count = std::ceil((hi[k] - lo[k] + 2 * margin) / spacing) + 1;
        if (!std::isfinite(count) || count > std::numeric_limits<int>::max())
            throw std::invalid_argument("grid too large");
        g.counts[k] = std::max(2, static_cast<int>(count));
    }
    g.steps = spacing * Eigen::Matrix3d::Identity();
    return g;
}

namespace {
void check_margin_and_spacing(double margin, double spacing) {
    if (!std::isfinite(margin) || margin < 0 || !std::isfinite(spacing) || spacing <= 0)
        throw std::invalid_argument("invalid grid margin or spacing");
}

// Fewest evenly spaced points, both ends included, covering `length` with a
// spacing that does not exceed `spacing`. A length that is a multiple of the
// spacing up to rounding (29.000000000000004 steps) keeps that multiple.
int points_covering(double length, double spacing) {
    const double steps = std::ceil(length / spacing * (1.0 - 1e-12));
    if (!std::isfinite(steps) || steps + 1 > std::numeric_limits<int>::max())
        throw std::invalid_argument("grid too large");
    return std::max(2, static_cast<int>(steps) + 1);
}
} // namespace

Grid Grid::line(const Eigen::Vector3d& start, const Eigen::Vector3d& end, double margin,
                double spacing) {
    check_margin_and_spacing(margin, spacing);
    if (!start.allFinite() || !end.allFinite())
        throw std::invalid_argument("invalid line: nonfinite end point");
    const double length = (end - start).norm();
    if (!(length > 0)) throw std::invalid_argument("invalid line: coincident end points");
    const Eigen::Vector3d direction = (end - start) / length;
    // Two unit vectors orthogonal to the line complete the grid axes; the
    // Cartesian axis least aligned with it keeps them well conditioned.
    Eigen::Index least = 0;
    direction.cwiseAbs().minCoeff(&least);
    const Eigen::Vector3d p = direction.cross(Eigen::Vector3d::Unit(least)).normalized();
    const double total = length + 2 * margin;
    Grid g;
    g.counts = Eigen::Vector3i(points_covering(total, spacing), 1, 1);
    g.origin = start - margin * direction;
    g.steps.col(0) = total / (g.counts[0] - 1) * direction;
    g.steps.col(1) = p;
    g.steps.col(2) = direction.cross(p);
    return g;
}

Grid Grid::plane(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                 double margin, double spacing) {
    check_margin_and_spacing(margin, spacing);
    if (!a.allFinite() || !b.allFinite() || !c.allFinite())
        throw std::invalid_argument("invalid plane: nonfinite point");
    const Eigen::Vector3d ab = b - a, ac = c - a;
    const double scale = ab.norm() * ac.norm();
    if (!(scale > 0) || ab.cross(ac).norm() <= 1e-8 * scale)
        throw std::invalid_argument("invalid plane: the three points are collinear");
    const Eigen::Vector3d u = ab.normalized();
    const Eigen::Vector3d v = (ac - ac.dot(u) * u).normalized();
    // In-plane coordinates of the three points: a (0, 0), b (|ab|, 0), c (cu, cv).
    const double cu = ac.dot(u), cv = ac.dot(v);
    const double u_low = std::min(0.0, cu) - margin, u_high = std::max(ab.norm(), cu) + margin;
    const double v_low = -margin, v_high = cv + margin;
    Grid g;
    g.counts = Eigen::Vector3i(points_covering(u_high - u_low, spacing),
                               points_covering(v_high - v_low, spacing), 1);
    g.origin = a + u_low * u + v_low * v;
    g.steps.col(0) = (u_high - u_low) / (g.counts[0] - 1) * u;
    g.steps.col(1) = (v_high - v_low) / (g.counts[1] - 1) * v;
    g.steps.col(2) = u.cross(v);
    return g;
}

Field::Field(Grid g, std::vector<Atom> atoms, std::vector<double> v)
    : grid_(std::move(g)), atoms_(std::move(atoms)), values_(std::move(v)) {
    if (values_.size() != detail::validated_point_count(grid_))
        throw std::invalid_argument("field size mismatch");
}

double Field::operator()(int i, int j, int k) const {
    if (i < 0 || j < 0 || k < 0 || i >= grid_.counts[0] || j >= grid_.counts[1] ||
        k >= grid_.counts[2])
        throw std::out_of_range("field index");
    return values_.at(detail::grid_index(grid_, i, j, k));
}
} // namespace chem::kumo
