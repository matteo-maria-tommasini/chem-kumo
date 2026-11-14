#include "chem/kumo/surface.hpp"
#include "electrostatic_potential.hpp"
#include "grid_geometry.hpp"
#include "isosurface.hpp"
#include "orbital_energy_fields.hpp"
#include "parallel_blocks.hpp"
#include "potential_clusters.hpp"

#include <Eigen/Geometry> // cross()

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace chem::kumo {

namespace detail {

/// Access to the Molecule internals that the surface functions need.
struct SurfaceAccess {
    /// Total density and gradient at the columns of `points`, with the
    /// default AO screening (a change of about 1e-13 e/bohr^3; about 1e-10
    /// bohr in the position of the rho = 0.001 isosurface).
    static std::vector<DensityDerivatives>
    total_density(const Molecule& molecule,
                  const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points) {
        return molecule
            .density_derivatives_shared(points, {Spin::Total}, 128, 1, default_screening_threshold)
            .at(Spin::Total);
    }

    static const BasisSet& basis(const Molecule& molecule) { return molecule.evaluator_.basis(); }

    static const BasisEvaluator& evaluator(const Molecule& molecule) { return molecule.evaluator_; }

    static const Eigen::MatrixXd& orbitals(const Molecule& molecule, Spin channel) {
        return molecule.mo_coefficients(channel);
    }
};

namespace {

constexpr double root_tolerance = 1e-10; // bohr, along the edge
constexpr int max_root_iterations = 200;

// The six tetrahedra of a cell around its main diagonal: corners 0,
// e_a, e_a + e_b and (1,1,1) for each permutation (a, b, c) of the axes.
// Corner masks have bit 0 for i, bit 1 for j, bit 2 for k.
constexpr std::array<std::array<int, 2>, 6> tetrahedron_axes = {
    {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}}};

} // namespace

bool isosurface_reaches_boundary(const Grid& grid, std::span<const double> samples,
                                 double isovalue) {
    const int nx = grid.counts[0], ny = grid.counts[1], nz = grid.counts[2];
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            const bool face = i == 0 || i == nx - 1 || j == 0 || j == ny - 1;
            for (int k = 0; k < nz; k += face || nz == 1 ? 1 : nz - 1)
                if (samples[grid_index(grid, i, j, k)] >= isovalue) return true;
        }
    return false;
}

Isosurface extract_isosurface(const Grid& grid, std::span<const double> samples, double isovalue,
                              const ScalarFunction& function) {
    const std::size_t npoints = validated_point_count(grid);
    if (grid.counts.minCoeff() < 2)
        throw std::invalid_argument("isosurface grid needs at least two points along each axis");
    if (samples.size() != npoints)
        throw std::invalid_argument("isosurface samples do not match the grid");
    if (!std::isfinite(isovalue)) throw std::invalid_argument("isovalue must be finite");
    if (!std::all_of(samples.begin(), samples.end(), [](double v) { return std::isfinite(v); }))
        throw std::invalid_argument("isosurface samples must be finite");
    if (isosurface_reaches_boundary(grid, samples, isovalue))
        throw std::invalid_argument(
            "the isosurface reaches the grid boundary, so it would not be closed");

    // One vertex per crossed grid edge, keyed by its lower end and direction.
    struct Edge {
        std::size_t inside = 0, outside = 0; // grid points (linear) on either side
        Eigen::Vector3i cell;                // lower end (i, j, k)
    };

    std::vector<Edge> edges;
    std::unordered_map<std::uint64_t, std::size_t> edge_ids;

    struct Polygon {
        std::array<std::size_t, 4> vertex{};
        int count = 0;
    };

    std::vector<Polygon> polygons;

    const auto inside = [&](std::size_t linear) { return samples[linear] >= isovalue; };
    for (int i = 0; i + 1 < grid.counts[0]; ++i)
        for (int j = 0; j + 1 < grid.counts[1]; ++j)
            for (int k = 0; k + 1 < grid.counts[2]; ++k) {
                std::array<std::size_t, 8> corner{};
                int ninside = 0;
                for (int mask = 0; mask < 8; ++mask) {
                    corner[static_cast<std::size_t>(mask)] =
                        grid_index(grid, i + (mask & 1), j + ((mask >> 1) & 1), k + (mask >> 2));
                    ninside += inside(corner[static_cast<std::size_t>(mask)]);
                }
                if (ninside == 0 || ninside == 8) continue;
                std::array<Eigen::Vector3d, 8> position;
                for (std::size_t mask = 0; mask < 8; ++mask)
                    position[mask] = grid_point(grid, corner[mask]);

                // Vertex on the edge from inside corner a to outside corner b.
                // The masks of a tetrahedron are nested, so the smaller one is
                // the lower end of the edge.
                const auto vertex = [&](int a, int b) {
                    const int low = std::min(a, b), high = std::max(a, b);
                    const auto base = corner[static_cast<std::size_t>(low)];
                    const std::uint64_t key = static_cast<std::uint64_t>(base) * 8 +
                                              static_cast<std::uint64_t>(high ^ low);
                    const auto [it, added] = edge_ids.try_emplace(key, edges.size());
                    if (added)
                        edges.push_back(
                            {corner[static_cast<std::size_t>(a)],
                             corner[static_cast<std::size_t>(b)],
                             Eigen::Vector3i(i + (low & 1), j + ((low >> 1) & 1), k + (low >> 2))});
                    return it->second;
                };

                for (const auto& axes : tetrahedron_axes) {
                    const int first = 1 << axes[0];
                    const std::array<int, 4> masks = {0, first, first | (1 << axes[1]), 7};
                    std::array<int, 4> in{}, out{};
                    int nin = 0, nout = 0;
                    for (int mask : masks) {
                        if (inside(corner[static_cast<std::size_t>(mask)]))
                            in[static_cast<std::size_t>(nin++)] = mask;
                        else
                            out[static_cast<std::size_t>(nout++)] = mask;
                    }
                    if (nin == 0 || nout == 0) continue;
                    // Polygon corners as (inside, outside) corner pairs, in cyclic order.
                    std::array<std::array<int, 2>, 4> ends{};
                    int count = 3;
                    if (nin == 1)
                        ends = {{{in[0], out[0]}, {in[0], out[1]}, {in[0], out[2]}, {}}};
                    else if (nout == 1)
                        ends = {{{in[0], out[0]}, {in[1], out[0]}, {in[2], out[0]}, {}}};
                    else {
                        ends = {
                            {{in[0], out[0]}, {in[0], out[1]}, {in[1], out[1]}, {in[1], out[0]}}};
                        count = 4;
                    }
                    // Orient with the edge midpoints, which depend only on the
                    // cell, so that the orientation is combinatorial: the
                    // normal points from the inside corners to the outside ones.
                    std::array<Eigen::Vector3d, 4> mid;
                    for (int c = 0; c < count; ++c) {
                        const auto& e = ends[static_cast<std::size_t>(c)];
                        mid[static_cast<std::size_t>(c)] =
                            0.5 * (position[static_cast<std::size_t>(e[0])] +
                                   position[static_cast<std::size_t>(e[1])]);
                    }
                    const Eigen::Vector3d normal =
                        count == 3 ? Eigen::Vector3d((mid[1] - mid[0]).cross(mid[2] - mid[0]))
                                   : Eigen::Vector3d((mid[2] - mid[0]).cross(mid[3] - mid[1]));
                    Eigen::Vector3d outward = Eigen::Vector3d::Zero();
                    for (int c = 0; c < nout; ++c)
                        outward +=
                            position[static_cast<std::size_t>(out[static_cast<std::size_t>(c)])] /
                            nout;
                    for (int c = 0; c < nin; ++c)
                        outward -=
                            position[static_cast<std::size_t>(in[static_cast<std::size_t>(c)])] /
                            nin;
                    if (normal.dot(outward) < 0) std::reverse(ends.begin(), ends.begin() + count);
                    Polygon polygon;
                    polygon.count = count;
                    for (int c = 0; c < count; ++c) {
                        const auto& e = ends[static_cast<std::size_t>(c)];
                        polygon.vertex[static_cast<std::size_t>(c)] = vertex(e[0], e[1]);
                    }
                    polygons.push_back(polygon);
                }
            }

    // Place every vertex on the isosurface: the root of f - isovalue along
    // its edge, x(t) = inside + t (outside - inside), bracketed by [0, 1].
    const std::size_t nvertices = edges.size();
    std::vector<double> t(nvertices), low(nvertices, 0.0), high(nvertices, 1.0);
    for (std::size_t v = 0; v < nvertices; ++v) {
        const double a = samples[edges[v].inside], b = samples[edges[v].outside];
        t[v] = std::clamp((a - isovalue) / (a - b), 0.0, 1.0);
    }
    const auto end_points = [&](std::size_t v) {
        return std::pair{grid_point(grid, edges[v].inside), grid_point(grid, edges[v].outside)};
    };
    std::vector<std::size_t> active(nvertices);
    std::iota(active.begin(), active.end(), std::size_t{0});
    std::vector<double> values;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points, gradients;
    for (int iteration = 0; !active.empty(); ++iteration) {
        if (iteration == max_root_iterations)
            throw std::runtime_error("isosurface vertices did not converge");
        const auto count = static_cast<Eigen::Index>(active.size());
        points.resize(3, count);
        gradients.setZero(3, count);
        values.assign(active.size(), 0.0);
        for (Eigen::Index p = 0; p < count; ++p) {
            const auto v = active[static_cast<std::size_t>(p)];
            const auto [a, b] = end_points(v);
            points.col(p) = a + t[v] * (b - a);
        }
        function(points, values, gradients);
        std::vector<std::size_t> next;
        for (Eigen::Index p = 0; p < count; ++p) {
            const auto v = active[static_cast<std::size_t>(p)];
            const double g = values[static_cast<std::size_t>(p)] - isovalue;
            if (!std::isfinite(g) || !gradients.col(p).allFinite())
                throw std::runtime_error("isosurface function is not finite");
            if (g == 0.0) continue;
            (g > 0.0 ? low[v] : high[v]) = t[v];
            const auto [a, b] = end_points(v);
            const double length = (b - a).norm();
            double step = t[v] - g / gradients.col(p).dot(b - a);
            if (!(step > low[v] && step < high[v])) step = 0.5 * (low[v] + high[v]);
            const double moved = std::abs(step - t[v]) * length;
            t[v] = step;
            if (moved > root_tolerance && (high[v] - low[v]) * length > root_tolerance)
                next.push_back(v);
        }
        active = std::move(next);
    }

    Isosurface surface;
    surface.mesh.vertices.resize(nvertices);
    surface.vertex_cells.resize(nvertices);
    for (std::size_t v = 0; v < nvertices; ++v) {
        const auto [a, b] = end_points(v);
        surface.mesh.vertices[v] = a + t[v] * (b - a);
        surface.vertex_cells[v] = edges[v].cell;
    }
    const auto& x = surface.mesh.vertices;
    for (const auto& polygon : polygons) {
        const auto& p = polygon.vertex;
        if (polygon.count == 3)
            surface.mesh.triangles.push_back({p[0], p[1], p[2]});
        else if ((x[p[0]] - x[p[2]]).squaredNorm() <= (x[p[1]] - x[p[3]]).squaredNorm()) {
            surface.mesh.triangles.push_back({p[0], p[1], p[2]});
            surface.mesh.triangles.push_back({p[0], p[2], p[3]});
        } else {
            surface.mesh.triangles.push_back({p[0], p[1], p[3]});
            surface.mesh.triangles.push_back({p[1], p[2], p[3]});
        }
    }
    return surface;
}

namespace {

// Total-density isosurface of density_isosurface(), with its vertex cells
// and the grid they refer to.
struct DensitySurface {
    Isosurface surface;
    Grid grid;
};

DensitySurface density_surface(const Molecule& molecule, double isovalue, double spacing) {
    if (!std::isfinite(isovalue) || isovalue < minimum_density_isovalue)
        throw std::invalid_argument("density isovalue must be finite and at least 1e-8 e/bohr^3");
    if (!std::isfinite(spacing) || spacing <= 0.0)
        throw std::invalid_argument("surface grid spacing must be positive and finite");
    const ScalarFunction density =
        [&molecule](const Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                    std::vector<double>& values,
                    Eigen::Matrix<double, 3, Eigen::Dynamic>& gradients) {
            const auto d = SurfaceAccess::total_density(molecule, points);
            for (std::size_t p = 0; p < d.size(); ++p) {
                values[p] = d[p].value;
                gradients.col(static_cast<Eigen::Index>(p)) = d[p].gradient;
            }
        };
    for (double margin = 4.0; margin <= 16.0; margin += 2.0) {
        const Grid grid = Grid::around(molecule.atoms(), margin, spacing);
        validated_point_count(grid);
        const auto field = molecule.densities_blocked(grid, {Spin::Total}).at(Spin::Total);
        if (isosurface_reaches_boundary(grid, field.values(), isovalue)) continue;
        return {extract_isosurface(grid, field.values(), isovalue, density), grid};
    }
    throw std::runtime_error("the density isosurface extends more than 16 bohr beyond the atoms; "
                             "use a larger isovalue");
}

// Calls piece(area, w0, w1, w2) for triangles covering the part of triangle x
// where the linear interpolation of its vertex values v is >= level (upper)
// or < level (lower); w are the interpolated values at their corners.
template <typename Piece>
void clip(const std::array<Eigen::Vector3d, 3>& x, const std::array<double, 3>& v, double level,
          bool upper, const Piece& piece) {
    std::array<Eigen::Vector3d, 4> px;
    std::array<double, 4> pv{};
    int n = 0;
    const auto keep = [&](double value) { return upper ? value >= level : value < level; };
    for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t j = (i + 1) % 3;
        const bool keep_i = keep(v[i]), keep_j = keep(v[j]);
        if (keep_i) {
            px[static_cast<std::size_t>(n)] = x[i];
            pv[static_cast<std::size_t>(n++)] = v[i];
        }
        if (keep_i != keep_j) {
            const double s = (level - v[i]) / (v[j] - v[i]);
            px[static_cast<std::size_t>(n)] = x[i] + s * (x[j] - x[i]);
            pv[static_cast<std::size_t>(n++)] = level;
        }
    }
    for (int i = 1; i + 1 < n; ++i) {
        const auto a = static_cast<std::size_t>(i), b = a + 1;
        piece(0.5 * (px[a] - px[0]).cross(px[b] - px[0]).norm(), pv[0], pv[a], pv[b]);
    }
}

// Integral over a triangle of area `area` of (w - s)^2, w linear with corner
// values w0, w1, w2.
double square_integral(double area, double w0, double w1, double w2, double s) {
    const double a = w0 - s, b = w1 - s, c = w2 - s;
    return area / 6.0 * (a * a + b * b + c * c + a * b + b * c + c * a);
}

} // namespace
} // namespace detail

double TriangleMesh::area() const {
    double total = 0.0;
    for (const auto& t : triangles) {
        const auto& a = vertices.at(t[0]);
        total += 0.5 * (vertices.at(t[1]) - a).cross(vertices.at(t[2]) - a).norm();
    }
    return total;
}

double TriangleMesh::volume() const {
    if (vertices.empty()) return 0.0;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (const auto& v : vertices)
        center += v;
    center /= static_cast<double>(vertices.size());
    double total = 0.0;
    for (const auto& t : triangles)
        total += (vertices.at(t[0]) - center)
                     .dot((vertices.at(t[1]) - center).cross(vertices.at(t[2]) - center));
    return total / 6.0;
}

SurfaceStatistics surface_statistics(const TriangleMesh& mesh, std::span<const double> values,
                                     double polar_threshold) {
    const std::size_t nvertices = mesh.vertices.size();
    if (values.size() != nvertices)
        throw std::invalid_argument("surface values do not match the mesh vertices");
    if (!std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }))
        throw std::invalid_argument("surface values must be finite");
    if (!std::isfinite(polar_threshold) || polar_threshold < 0.0)
        throw std::invalid_argument("polar threshold must be finite and nonnegative");
    for (const auto& t : mesh.triangles)
        if (t[0] >= nvertices || t[1] >= nvertices || t[2] >= nvertices)
            throw std::invalid_argument("surface triangle refers to a missing vertex");

    SurfaceStatistics s;
    const auto corners = [&](const std::array<std::size_t, 3>& t) {
        return std::pair{std::array<Eigen::Vector3d, 3>{mesh.vertices[t[0]], mesh.vertices[t[1]],
                                                        mesh.vertices[t[2]]},
                         std::array<double, 3>{values[t[0]], values[t[1]], values[t[2]]}};
    };
    const auto mean = [](double a, double b, double c) { return (a + b + c) / 3.0; };

    double integral = 0.0, positive_integral = 0.0, negative_integral = 0.0;
    double above_area = 0.0, below_area = 0.0; // V > threshold, V < -threshold
    for (const auto& t : mesh.triangles) {
        const auto [x, v] = corners(t);
        detail::clip(x, v, 0.0, true, [&](double a, double w0, double w1, double w2) {
            s.positive_area += a;
            positive_integral += a * mean(w0, w1, w2);
        });
        detail::clip(x, v, 0.0, false, [&](double a, double w0, double w1, double w2) {
            s.negative_area += a;
            negative_integral += a * mean(w0, w1, w2);
        });
        detail::clip(x, v, polar_threshold, true,
                     [&](double a, double, double, double) { above_area += a; });
        detail::clip(x, v, -polar_threshold, false,
                     [&](double a, double, double, double) { below_area += a; });
    }
    s.area = s.positive_area + s.negative_area;
    integral = positive_integral + negative_integral;
    if (!(s.area > 0.0) || !std::isfinite(s.area))
        throw std::invalid_argument("surface area must be positive");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    s.average = integral / s.area;
    s.positive_average = s.positive_area > 0.0 ? positive_integral / s.positive_area : nan;
    s.negative_average = s.negative_area > 0.0 ? negative_integral / s.negative_area : nan;
    s.polarity_index = (positive_integral - negative_integral) / s.area;
    s.nonpolar_area = std::max(0.0, s.area - above_area - below_area);

    double positive_square = 0.0, negative_square = 0.0, deviation = 0.0;
    for (const auto& t : mesh.triangles) {
        const auto [x, v] = corners(t);
        if (s.positive_area > 0.0)
            detail::clip(x, v, 0.0, true, [&](double a, double w0, double w1, double w2) {
                positive_square += detail::square_integral(a, w0, w1, w2, s.positive_average);
            });
        if (s.negative_area > 0.0)
            detail::clip(x, v, 0.0, false, [&](double a, double w0, double w1, double w2) {
                negative_square += detail::square_integral(a, w0, w1, w2, s.negative_average);
            });
        detail::clip(x, v, s.average, true, [&](double a, double w0, double w1, double w2) {
            deviation += a * (mean(w0, w1, w2) - s.average);
        });
        detail::clip(x, v, s.average, false, [&](double a, double w0, double w1, double w2) {
            deviation += a * (s.average - mean(w0, w1, w2));
        });
    }
    s.positive_variance = s.positive_area > 0.0 ? positive_square / s.positive_area : 0.0;
    s.negative_variance = s.negative_area > 0.0 ? negative_square / s.negative_area : 0.0;
    s.total_variance = s.positive_variance + s.negative_variance;
    s.balance = s.total_variance > 0.0 ? s.positive_variance * s.negative_variance /
                                             (s.total_variance * s.total_variance)
                                       : nan;
    s.charge_separation = deviation / s.area;

    const auto [lowest, highest] = std::minmax_element(values.begin(), values.end());
    s.minimum_vertex = static_cast<std::size_t>(lowest - values.begin());
    s.maximum_vertex = static_cast<std::size_t>(highest - values.begin());
    s.minimum = *lowest;
    s.maximum = *highest;

    // Local extrema: below (above) every neighbour and neighbour of a neighbour.
    std::vector<std::vector<std::size_t>> neighbours(nvertices);
    for (const auto& t : mesh.triangles)
        for (std::size_t a = 0; a < 3; ++a) {
            neighbours[t[a]].push_back(t[(a + 1) % 3]);
            neighbours[t[a]].push_back(t[(a + 2) % 3]);
        }
    for (auto& list : neighbours) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
    const auto extremum = [&](std::size_t i, auto beyond) {
        if (neighbours[i].empty()) return false;
        for (std::size_t j : neighbours[i]) {
            if (j == i) continue;
            if (!beyond(values[i], values[j])) return false;
            for (std::size_t k : neighbours[j])
                if (k != i && !beyond(values[i], values[k])) return false;
        }
        return true;
    };
    for (std::size_t i = 0; i < nvertices; ++i) {
        if (extremum(i, std::less<double>{})) s.local_minima.push_back(i);
        if (extremum(i, std::greater<double>{})) s.local_maxima.push_back(i);
    }
    std::sort(s.local_minima.begin(), s.local_minima.end(), [&](std::size_t a, std::size_t b) {
        return values[a] < values[b] || (values[a] == values[b] && a < b);
    });
    std::sort(s.local_maxima.begin(), s.local_maxima.end(), [&](std::size_t a, std::size_t b) {
        return values[a] > values[b] || (values[a] == values[b] && a < b);
    });
    return s;
}

TriangleMesh density_isosurface(const Molecule& molecule, double isovalue, double spacing) {
    return detail::density_surface(molecule, isovalue, spacing).surface.mesh;
}

std::vector<double> average_local_ionization_energy(
    const Molecule& molecule,
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points) {
    if (!points.allFinite()) throw std::invalid_argument("ionization energy points must be finite");
    const detail::OrbitalEnergyFields fields(molecule, detail::SurfaceAccess::evaluator(molecule),
                                             detail::SurfaceAccess::orbitals(molecule, Spin::Alpha),
                                             detail::SurfaceAccess::orbitals(molecule, Spin::Beta),
                                             false, "the average local ionization energy");
    const auto count = static_cast<std::size_t>(points.cols());
    std::vector<double> result(count);
    constexpr std::size_t block_size = 128;
    detail::for_each_point_block(
        count, block_size, "too many ionization energy blocks for parallel evaluation",
        [&](detail::PointBlock range) {
            const std::size_t start = range.start;
            const auto b = static_cast<Eigen::Index>(range.size);
            const auto sums =
                fields.evaluate(points.middleCols(static_cast<Eigen::Index>(start), b));
            for (std::size_t t = 0; t < static_cast<std::size_t>(b); ++t) {
                const double density = sums.occupied_density[t];
                if (!(density > 0.0))
                    throw std::invalid_argument(
                        "the average local ionization energy is undefined where the density "
                        "vanishes");
                result[start + t] = sums.occupied_energy[t] / density;
            }
        });
    return result;
}

IonizationEnergySurface local_ionization_energy_surface(const Molecule& molecule, double isovalue,
                                                        double spacing) {
    // Check the orbitals before building the surface (no points to evaluate).
    (void)average_local_ionization_energy(molecule, Eigen::Matrix<double, 3, 0>());
    auto surface = detail::density_surface(molecule, isovalue, spacing).surface;
    auto& mesh = surface.mesh;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(
        3, static_cast<Eigen::Index>(mesh.vertices.size()));
    for (std::size_t v = 0; v < mesh.vertices.size(); ++v)
        points.col(static_cast<Eigen::Index>(v)) = mesh.vertices[v];
    IonizationEnergySurface result;
    result.ionization_energy = average_local_ionization_energy(molecule, points);
    result.statistics = surface_statistics(mesh, result.ionization_energy);
    result.mesh = std::move(mesh);
    return result;
}

PotentialSurface electrostatic_potential_surface(const Molecule& molecule,
                                                 const SurfaceOptions& options) {
    if (!std::isfinite(options.potential_tolerance) || options.potential_tolerance < 0.0)
        throw std::invalid_argument("potential tolerance must be finite and nonnegative");
    auto density = detail::density_surface(molecule, options.isovalue, options.spacing);
    auto& surface = density.surface;
    const Grid& grid = density.grid;
    const detail::ElectrostaticPotential potential(detail::SurfaceAccess::basis(molecule),
                                                   molecule.density_matrix(Spin::Total),
                                                   molecule.atoms());

    // Vertices are evaluated in clusters: those whose edges start in the same
    // block of tile^3 grid cells, so they lie within the index box
    // [tile T, tile (T + 1)] around the block centre tile (T + 1/2).
    constexpr int tile = 8;
    const auto& mesh = surface.mesh;
    const std::size_t nvertices = mesh.vertices.size();
    const auto tile_of = [&](std::size_t v) {
        const Eigen::Vector3i& c = surface.vertex_cells[v];
        return std::array<int, 3>{c[0] / tile, c[1] / tile, c[2] / tile};
    };
    const auto clusters = detail::cluster_by(nvertices, tile_of);
    const auto plan = detail::cluster_plan(
        potential, detail::grid_box_radius(grid.steps, Eigen::Vector3d::Constant(tile)),
        options.potential_tolerance);
    const auto gather = [&](std::size_t c, Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                            std::vector<std::size_t>& targets) {
        detail::gather_cluster(
            clusters, c, [&](std::size_t v) { return mesh.vertices[v]; }, points, targets);
        const auto block = tile_of(targets.front());
        return Eigen::Vector3d(
            grid.origin +
            grid.steps *
                (tile * (Eigen::Vector3d(block[0], block[1], block[2]).array() + 0.5)).matrix());
    };

    PotentialSurface result;
    result.potential.assign(nvertices, 0.0);
    detail::evaluate_in_clusters(potential, plan, clusters.size(), gather, result.potential,
                                 nvertices, {});
    result.statistics = surface_statistics(mesh, result.potential);
    result.mesh = std::move(surface.mesh);
    return result;
}

MappedSurface map_on_density_isosurface(const Molecule& molecule,
                                        const std::vector<Property>& properties, double isovalue,
                                        double spacing, const EvaluationOptions& options) {
    if (properties.empty()) throw std::invalid_argument("no property to map on the surface");
    MappedSurface result;
    Request request;
    for (const Property& property : properties) {
        if (property.kind == PropertyKind::Gradient || property.kind == PropertyKind::Hessian)
            throw std::invalid_argument(
                "only one-component properties can be mapped on a surface (not the gradient "
                "or the Hessian)");
        if (std::find(result.properties.begin(), result.properties.end(), property) ==
            result.properties.end()) {
            result.properties.push_back(property);
            request.add(property);
        }
    }
    // Check the request (components, orbitals) before building the surface.
    const Evaluator evaluator(molecule, options);
    (void)evaluator.evaluate(Eigen::Matrix<double, 3, Eigen::Dynamic>(3, 0), request);
    result.mesh = detail::density_surface(molecule, isovalue, spacing).surface.mesh;
    const auto& vertices = result.mesh.vertices;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(vertices.size()));
    for (std::size_t v = 0; v < vertices.size(); ++v)
        points.col(static_cast<Eigen::Index>(v)) = vertices[v];
    const Results results = evaluator.evaluate(points, request);
    result.backend_used = results.backend_used();
    for (const Property& property : result.properties) {
        const auto values = results.get(property).component();
        result.values.emplace_back(values.begin(), values.end());
        result.statistics.push_back(surface_statistics(result.mesh, result.values.back(), 0.0));
    }
    return result;
}

void write_ply(const std::filesystem::path& path, const TriangleMesh& mesh,
               const std::vector<std::pair<std::string, std::vector<double>>>& fields,
               const std::vector<std::string>& comments) {
    std::vector<std::string> names{"x", "y", "z"};
    for (const auto& [name, values] : fields) {
        if (name.empty() || !std::all_of(name.begin(), name.end(), [](unsigned char c) {
                return std::isalnum(c) || c == '_';
            }))
            throw std::invalid_argument("invalid PLY property name '" + name + "'");
        if (std::find(names.begin(), names.end(), name) != names.end())
            throw std::invalid_argument("duplicate PLY property name '" + name + "'");
        names.push_back(name);
        if (values.size() != mesh.vertices.size())
            throw std::invalid_argument("PLY property '" + name + "' needs one value per vertex");
        if (!std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }))
            throw std::invalid_argument("nonfinite value in PLY property '" + name + "'");
    }
    for (const auto& comment : comments)
        if (comment.find_first_of("\r\n") != std::string::npos)
            throw std::invalid_argument("PLY comments must not contain line breaks");
    for (const auto& v : mesh.vertices)
        if (!v.allFinite()) throw std::invalid_argument("nonfinite PLY vertex");
    for (const auto& t : mesh.triangles)
        for (const std::size_t v : t)
            if (v >= mesh.vertices.size())
                throw std::invalid_argument("PLY triangle refers to a missing vertex");

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open PLY output " + path.string());
    out << "ply\nformat ascii 1.0\n";
    for (const auto& comment : comments)
        out << "comment " << comment << '\n';
    out << "element vertex " << mesh.vertices.size() << '\n';
    for (const auto& name : names)
        out << "property double " << name << '\n';
    out << "element face " << mesh.triangles.size() << '\n'
        << "property list uchar int vertex_indices\nend_header\n";
    out << std::scientific << std::setprecision(9);
    for (std::size_t v = 0; v < mesh.vertices.size(); ++v) {
        out << mesh.vertices[v].x() << ' ' << mesh.vertices[v].y() << ' ' << mesh.vertices[v].z();
        for (const auto& field : fields)
            out << ' ' << field.second[v];
        out << '\n';
    }
    for (const auto& t : mesh.triangles)
        out << "3 " << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
    out.flush();
    if (!out) throw std::runtime_error("PLY write failed: " + path.string());
}

} // namespace chem::kumo
