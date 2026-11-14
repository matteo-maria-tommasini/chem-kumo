#pragma once

// Electrostatic potential evaluated in clusters of nearby points. Grids use
// tiles of grid points, isosurfaces tiles of the cells that hold the
// vertices, and arbitrary points fixed cubes in space; in every case each
// cluster is evaluated around its own centre with one far-field plan shared by
// all clusters. Clusters depend only on the points, and each value only on its
// point, its cluster centre and the plan, so the result does not depend on the
// number of threads.

#include "chem/kumo/density.hpp"
#include "electrostatic_potential.hpp"
#include "parallel_blocks.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

namespace chem::kumo::detail {

/// Far-field plan for clusters whose points lie within `radius` of their
/// centre, slightly enlarged to cover rounding in the point coordinates.
inline ElectrostaticPotential::FarFieldPlan cluster_plan(const ElectrostaticPotential& potential,
                                                         double radius, double tolerance) {
    return potential.far_field_plan(radius * (1.0 + 1e-10) + 1e-12, tolerance);
}

/// Largest distance from the centre of a box of `edges[a]` grid steps along
/// each grid axis to its corners.
inline double grid_box_radius(const Eigen::Matrix3d& steps, const Eigen::Vector3d& edges) {
    double radius = 0.0;
    for (int corner = 0; corner < 8; ++corner) {
        Eigen::Vector3d offset;
        for (int a = 0; a < 3; ++a)
            offset[a] = ((corner >> a) & 1 ? 0.5 : -0.5) * edges[a];
        radius = std::max(radius, (steps * offset).norm());
    }
    return radius;
}

/// Points 0, ..., n - 1 grouped by a key: cluster c holds the points
/// order[starts[c]], ..., order[starts[c + 1] - 1], in increasing index.
struct PointClusters {
    std::vector<std::size_t> order;
    std::vector<std::size_t> starts;

    [[nodiscard]] std::size_t size() const { return starts.size() - 1; }
};

/// Group points 0, ..., n - 1 by key(p) (any type with operator<), clusters
/// in increasing key order.
template <typename Key> PointClusters cluster_by(std::size_t n, const Key& key) {
    PointClusters clusters;
    clusters.order.resize(n);
    std::iota(clusters.order.begin(), clusters.order.end(), std::size_t{0});
    std::stable_sort(clusters.order.begin(), clusters.order.end(),
                     [&](std::size_t a, std::size_t b) { return key(a) < key(b); });
    for (std::size_t p = 0; p < n; ++p)
        if (p == 0 || key(clusters.order[p]) != key(clusters.order[p - 1]))
            clusters.starts.push_back(p);
    clusters.starts.push_back(n);
    return clusters;
}

/// Potential at `total` points in `nclusters` clusters, on OpenMP threads.
/// gather(c, points, targets) fills the points of cluster c (columns) and
/// their indices in `output`, and returns the cluster centre; the points must
/// lie within plan.radius of it.
template <typename Gather>
void evaluate_in_clusters(const ElectrostaticPotential& potential,
                          const ElectrostaticPotential::FarFieldPlan& plan, std::size_t nclusters,
                          const Gather& gather, std::vector<double>& output, std::size_t total,
                          const std::function<void(std::size_t, std::size_t)>& progress) {
    Progress report(progress, total);
    for_each_block(nclusters, [&](std::size_t c) {
        Eigen::Matrix<double, 3, Eigen::Dynamic> points;
        std::vector<std::size_t> targets;
        const Eigen::Vector3d center = gather(c, points, targets);
        std::vector<double> values(targets.size());
        potential.evaluate_cluster(points, center, plan, values);
        for (std::size_t p = 0; p < values.size(); ++p)
            output[targets[p]] = values[p];
        report.advance(values.size());
    });
}

/// gather() helper for PointClusters: copies the points of cluster c, with
/// point_at(p) the coordinates of point p, and their indices.
template <typename PointAt>
void gather_cluster(const PointClusters& clusters, std::size_t c, const PointAt& point_at,
                    Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                    std::vector<std::size_t>& targets) {
    const std::size_t begin = clusters.starts[c], end = clusters.starts[c + 1];
    targets.assign(clusters.order.begin() + static_cast<std::ptrdiff_t>(begin),
                   clusters.order.begin() + static_cast<std::ptrdiff_t>(end));
    points.resize(3, static_cast<Eigen::Index>(targets.size()));
    for (std::size_t p = 0; p < targets.size(); ++p)
        points.col(static_cast<Eigen::Index>(p)) = point_at(targets[p]);
}

} // namespace chem::kumo::detail
