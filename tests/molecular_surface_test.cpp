// Molecular surfaces: the isosurface construction on analytic functions, the
// surface statistics on a mesh with exact answers, and the density
// isosurface and electrostatic potential surface of a molecule.

#include "chem/kumo.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "electrostatic_potential.hpp"
#include "isosurface.hpp"
#include "omp_support.hpp"
#include "test_support.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using chem::kumo::SurfaceStatistics;
using chem::kumo::TriangleMesh;
using Points = Eigen::Matrix<double, 3, Eigen::Dynamic>;

// Every directed edge appears once and its reverse once (closed, consistently
// oriented 2-manifold); returns the Euler characteristic V - E + F.
long check_closed(const TriangleMesh& mesh, const std::string& label) {
    std::map<std::pair<std::size_t, std::size_t>, int> edges;
    for (const auto& t : mesh.triangles)
        for (std::size_t a = 0; a < 3; ++a)
            ++edges[{t[a], t[(a + 1) % 3]}];
    std::size_t bad = 0;
    for (const auto& [edge, count] : edges)
        if (count != 1 || !edges.contains({edge.second, edge.first})) ++bad;
    kumo_test::check(bad == 0, label +
                                   ": every edge shared by two triangles with opposite "
                                   "orientations (" +
                                   std::to_string(bad) + " directed edges are not)");
    return static_cast<long>(mesh.vertices.size()) - static_cast<long>(edges.size() / 2) +
           static_cast<long>(mesh.triangles.size());
}

// Axis-aligned grid with the given spacing covering [-half, half]^3, with an
// offset so that no grid point is at the centre.
chem::kumo::Grid cube_grid(double half, double spacing) {
    chem::kumo::Grid grid;
    grid.origin = Eigen::Vector3d::Constant(-half) + Eigen::Vector3d(0.013, 0.007, 0.021);
    const int n = static_cast<int>(std::ceil(2 * half / spacing)) + 1;
    grid.counts = Eigen::Vector3i::Constant(n);
    grid.steps = spacing * Eigen::Matrix3d::Identity();
    return grid;
}

// Sphere |r| = radius as the isosurface f = -|r| = -radius.
TriangleMesh sphere(double radius, double spacing) {
    const auto grid = cube_grid(radius + 1.0, spacing);
    std::vector<double> samples;
    for (int i = 0; i < grid.counts[0]; ++i)
        for (int j = 0; j < grid.counts[1]; ++j)
            for (int k = 0; k < grid.counts[2]; ++k)
                samples.push_back(-(grid.origin + grid.steps * Eigen::Vector3d(i, j, k)).norm());
    const chem::kumo::detail::ScalarFunction f =
        [](const Points& points, std::vector<double>& values, Points& gradients) {
            for (Eigen::Index p = 0; p < points.cols(); ++p) {
                const double r = points.col(p).norm();
                values[static_cast<std::size_t>(p)] = -r;
                gradients.col(p) = -points.col(p) / r;
            }
        };
    return chem::kumo::detail::extract_isosurface(grid, samples, -radius, f).mesh;
}

void test_sphere() {
    const double radius = 1.7;
    const Eigen::Vector3d direction = Eigen::Vector3d(0.3, 0.2, 0.9).normalized();
    double previous_area_error = 0.0, previous_volume_error = 0.0;
    for (double spacing : {0.2, 0.1}) {
        const std::string label = "sphere, spacing " + kumo_test::num(spacing);
        const auto mesh = sphere(radius, spacing);
        kumo_test::check(check_closed(mesh, label) == 2, label + ": Euler characteristic 2");
        double off = 0.0;
        for (const auto& v : mesh.vertices)
            off = std::max(off, std::abs(v.norm() - radius));
        kumo_test::check(off <= 1e-9, label +
                                          ": vertices on the sphere within 1e-9 bohr (largest "
                                          "deviation " +
                                          kumo_test::num(off) + ")");
        // Inscribed polyhedron: area and volume converge to the sphere's from
        // below, at second order in the spacing.
        const double area = mesh.area(), volume = mesh.volume();
        const double exact_area = 4 * std::numbers::pi * radius * radius;
        const double exact_volume = 4.0 / 3.0 * std::numbers::pi * radius * radius * radius;
        const double area_error = (exact_area - area) / exact_area;
        const double volume_error = (exact_volume - volume) / exact_volume;
        kumo_test::check(area_error > 0 && area_error < 3e-3,
                         label + ": area below 4 pi R^2 by less than 0.3% (relative error " +
                             kumo_test::num(area_error) + ")");
        kumo_test::check(volume_error > 0 && volume_error < 6e-3,
                         label + ": volume below 4/3 pi R^3 by less than 0.6% (relative error " +
                             kumo_test::num(volume_error) + ")");
        if (previous_area_error > 0)
            kumo_test::check(previous_area_error / area_error > 3 &&
                                 previous_volume_error / volume_error > 3,
                             label +
                                 ": halving the spacing reduces the area and volume errors "
                                 "more than threefold (ratios " +
                                 kumo_test::num(previous_area_error / area_error) + ", " +
                                 kumo_test::num(previous_volume_error / volume_error) + ")");
        previous_area_error = area_error;
        previous_volume_error = volume_error;

        // V = n.r on a sphere: averages over hemispheres R/2, variance
        // R^2/3 - R^2/4 = R^2/12, nu = 1/4, Pi = MPI = R/2, and the band
        // |V| <= c has area 4 pi R c (Archimedes).
        std::vector<double> values;
        for (const auto& v : mesh.vertices)
            values.push_back(direction.dot(v));
        const double threshold = 0.4;
        const auto s = chem::kumo::surface_statistics(mesh, values, threshold);
        // Discretization errors at spacing 0.2: at most 1.4e-3 (relative).
        const double rel = 4e-3;
        kumo_test::check_close(s.area, area, 1e-12, 0, label + ": statistics area = mesh area");
        kumo_test::check_close(s.positive_area / s.area, 0.5, rel, rel, label + ": A+ / A");
        kumo_test::check(std::abs(s.average) < 2e-4 * radius,
                         label + ": average of n.r near 0 (got " + kumo_test::num(s.average) + ")");
        kumo_test::check_close(s.positive_average, radius / 2, rel, 0,
                               label + ": positive average");
        kumo_test::check_close(s.negative_average, -radius / 2, rel, 0,
                               label + ": negative average");
        kumo_test::check_close(s.positive_variance, radius * radius / 12, rel, 0,
                               label + ": positive variance");
        kumo_test::check_close(s.negative_variance, radius * radius / 12, rel, 0,
                               label + ": negative variance");
        kumo_test::check_close(s.balance, 0.25, 0, 1e-6, label + ": balance nu");
        kumo_test::check_close(s.charge_separation, radius / 2, rel, 0, label + ": Pi");
        kumo_test::check_close(s.polarity_index, radius / 2, rel, 0, label + ": MPI");
        kumo_test::check_close(s.nonpolar_area, 4 * std::numbers::pi * radius * threshold, rel, 0,
                               label + ": nonpolar area");
        // Local extrema: one minimum at -R n and one maximum at R n, with the
        // global extrema.
        if (kumo_test::check(s.local_minima.size() == 1 && s.local_maxima.size() == 1,
                             label + ": one local minimum and one local maximum (got " +
                                 std::to_string(s.local_minima.size()) + " and " +
                                 std::to_string(s.local_maxima.size()) + ")")) {
            kumo_test::check(s.local_minima[0] == s.minimum_vertex &&
                                 s.local_maxima[0] == s.maximum_vertex,
                             label + ": local extrema are the global ones");
            const double distance = (mesh.vertices[s.maximum_vertex] - radius * direction).norm();
            kumo_test::check(distance < spacing, label + ": maximum within one spacing of R n (" +
                                                     kumo_test::num(distance) + " bohr)");
            kumo_test::check_close(s.maximum, radius, 0, spacing * spacing, label + ": maximum");
            kumo_test::check_close(s.minimum, -radius, 0, spacing * spacing, label + ": minimum");
        }
    }
}

// Regular octahedron with vertices +-e_x, +-e_y, +-e_z and V = x. Faces have
// area sqrt(3)/2; the hand-computed integrals of the linear interpolation are
// exact, so the statistics agree to rounding.
TriangleMesh octahedron() {
    TriangleMesh mesh;
    mesh.vertices = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (std::size_t x : {0u, 1u})
        for (std::size_t y : {2u, 3u})
            for (std::size_t z : {4u, 5u}) {
                std::array<std::size_t, 3> t = {x, y, z};
                const auto& v = mesh.vertices;
                if ((v[y] - v[x]).cross(v[z] - v[x]).dot(v[x]) < 0) std::swap(t[1], t[2]);
                mesh.triangles.push_back(t);
            }
    return mesh;
}

void test_octahedron() {
    const auto mesh = octahedron();
    const double r3 = std::sqrt(3.0);
    kumo_test::check(check_closed(mesh, "octahedron") == 2, "octahedron: Euler characteristic 2");
    kumo_test::check_close(mesh.area(), 4 * r3, 1e-14, 0, "octahedron area");
    kumo_test::check_close(mesh.volume(), 4.0 / 3.0, 1e-14, 0, "octahedron volume");

    std::vector<double> x;
    for (const auto& v : mesh.vertices)
        x.push_back(v.x());
    // Positive half: 4 faces with centroid x = 1/3, so A+ = 2 sqrt(3) and the
    // average 1/3; integral of x^2 over a face (sqrt(3)/2)/6, so the mean of x^2
    // is 1/6 and the variance 1/6 - 1/9 = 1/18. |x| > 1/2 cuts a corner of area
    // (1/2)^2 sqrt(3)/2 from each face.
    const auto s = chem::kumo::surface_statistics(mesh, x, 0.5);
    const double t = 1e-14;
    kumo_test::check_close(s.area, 4 * r3, t, 0, "octahedron A");
    kumo_test::check_close(s.positive_area, 2 * r3, t, 0, "octahedron A+");
    kumo_test::check_close(s.negative_area, 2 * r3, t, 0, "octahedron A-");
    kumo_test::check_close(s.average, 0, 0, t, "octahedron average");
    kumo_test::check_close(s.positive_average, 1.0 / 3, t, 0, "octahedron positive average");
    kumo_test::check_close(s.negative_average, -1.0 / 3, t, 0, "octahedron negative average");
    kumo_test::check_close(s.positive_variance, 1.0 / 18, t, 0, "octahedron positive variance");
    kumo_test::check_close(s.negative_variance, 1.0 / 18, t, 0, "octahedron negative variance");
    kumo_test::check_close(s.total_variance, 1.0 / 9, t, 0, "octahedron total variance");
    kumo_test::check_close(s.balance, 0.25, t, 0, "octahedron nu");
    kumo_test::check_close(s.charge_separation, 1.0 / 3, t, 0, "octahedron Pi");
    kumo_test::check_close(s.polarity_index, 1.0 / 3, t, 0, "octahedron MPI");
    kumo_test::check_close(s.nonpolar_area, 3 * r3, t, 0, "octahedron nonpolar area");
    kumo_test::check(s.minimum_vertex == 1 && s.maximum_vertex == 0 && s.minimum == -1 &&
                         s.maximum == 1,
                     "octahedron: global minimum at -e_x and maximum at +e_x");
    kumo_test::check(s.local_minima == std::vector<std::size_t>{1} &&
                         s.local_maxima == std::vector<std::size_t>{0},
                     "octahedron: the only local extrema are -e_x and +e_x");

    // Shifted values: V = x + 1/4 crosses zero at x = -1/4. The positive
    // part is everything but the corners x < -1/4 (scale 3/4 of the faces at
    // -e_x), so A- = 4 (3/4)^2 sqrt(3)/2 = 9 sqrt(3) / 8.
    std::vector<double> shifted;
    for (double v : x)
        shifted.push_back(v + 0.25);
    const auto u = chem::kumo::surface_statistics(mesh, shifted);
    kumo_test::check_close(u.negative_area, 9 * r3 / 8, t, 0, "shifted octahedron A-");
    kumo_test::check_close(u.average, 0.25, t, 0, "shifted octahedron average");
    // On a corner triangle (apex value -3/4, base value 0) the mean is -1/4.
    kumo_test::check_close(u.negative_average, -0.25, t, 0, "shifted octahedron negative average");

    // All positive: no negative region; constant: no variance.
    std::vector<double> positive, constant(x.size(), 2.0);
    for (double v : x)
        positive.push_back(v + 2.0);
    const auto p = chem::kumo::surface_statistics(mesh, positive);
    kumo_test::check(p.negative_area == 0 && std::isnan(p.negative_average) &&
                         p.negative_variance == 0 && p.balance == 0,
                     "all-positive values: A- = 0, negative average NaN, variance 0, nu 0");
    const auto c = chem::kumo::surface_statistics(mesh, constant);
    kumo_test::check(c.total_variance == 0 && std::isnan(c.balance) && c.nonpolar_area == 0 &&
                         c.local_minima.empty() && c.local_maxima.empty(),
                     "constant values: variance 0, nu NaN, no nonpolar area, no local extrema");

    // Invalid input.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)chem::kumo::surface_statistics(mesh, std::vector<double>(5, 0.0)); },
        "surface_statistics rejects a value count different from the vertex count");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto bad = x;
            bad[2] = nan;
            (void)chem::kumo::surface_statistics(mesh, bad);
        },
        "surface_statistics rejects a NaN value");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)chem::kumo::surface_statistics(mesh, x, -1.0); },
        "surface_statistics rejects a negative polar threshold");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto bad = mesh;
            bad.triangles[0][1] = 6;
            (void)chem::kumo::surface_statistics(bad, x);
        },
        "surface_statistics rejects a triangle with a missing vertex");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto flat = mesh;
            for (auto& v : flat.vertices)
                v = Eigen::Vector3d::Zero();
            (void)chem::kumo::surface_statistics(flat, x);
        },
        "surface_statistics rejects a mesh of zero area");
}

void test_extraction_errors() {
    const auto grid = cube_grid(2.0, 0.5);
    const std::size_t n = static_cast<std::size_t>(grid.counts.prod());
    const chem::kumo::detail::ScalarFunction f = [](const Points&, std::vector<double>&, Points&) {
    };
    std::vector<double> samples(n, -1.0);
    samples[n / 2] = 1.0;
    using chem::kumo::detail::extract_isosurface;
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)extract_isosurface(grid, std::vector<double>(n - 1, -1.0), 0.0, f); },
        "extract_isosurface rejects a sample count different from the grid");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto bad = samples;
            bad[3] = std::numeric_limits<double>::infinity();
            (void)extract_isosurface(grid, bad, 0.0, f);
        },
        "extract_isosurface rejects a nonfinite sample");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)extract_isosurface(grid, samples, std::nan(""), f); },
        "extract_isosurface rejects a NaN isovalue");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto bad = samples;
            bad.back() = 1.0; // a corner of the grid
            (void)extract_isosurface(grid, bad, 0.0, f);
        },
        "extract_isosurface rejects a region that reaches the grid boundary");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto flat = grid;
            flat.counts[2] = 1;
            (void)extract_isosurface(flat, std::vector<double>(n / grid.counts[2], -1.0), 0.0, f);
        },
        "extract_isosurface rejects a grid with one point along an axis");
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            auto bad = grid;
            bad.steps.col(1) = bad.steps.col(0); // linearly dependent axes
            (void)extract_isosurface(bad, samples, 0.0, f);
        },
        "extract_isosurface rejects an invalid grid (detail::validated_point_count)");
    // Nothing above the isovalue: an empty mesh.
    const auto empty = extract_isosurface(grid, std::vector<double>(n, -1.0), 0.0, f);
    kumo_test::check(empty.mesh.vertices.empty() && empty.mesh.triangles.empty(),
                     "extract_isosurface gives an empty mesh when no sample reaches the isovalue");
}

void test_molecule(const char* path) {
    const std::string name = path;
    const auto mol = chem::kumo::Molecule::from_fchk(path);
    kumo_test::set_threads(1);
    const auto reference = chem::kumo::electrostatic_potential_surface(mol);
    const auto& mesh = reference.mesh;
    check_closed(mesh, name);
    kumo_test::check(mesh.volume() > 0, name + ": positive enclosed volume");

    // Vertices on the density isosurface, and the same mesh as density_isosurface().
    Points points(3, static_cast<Eigen::Index>(mesh.vertices.size()));
    for (std::size_t v = 0; v < mesh.vertices.size(); ++v)
        points.col(static_cast<Eigen::Index>(v)) = mesh.vertices[v];
    const auto density = mol.density_derivatives_blocked(points, chem::kumo::Spin::Total, 128, 0.0);
    double worst = 0.0;
    for (const auto& d : density)
        worst = std::max(worst, std::abs(d.value - 1e-3) / d.gradient.norm());
    kumo_test::check(worst <= 1e-9, name +
                                        ": vertices within 1e-9 bohr of rho = 0.001 (largest "
                                        "distance estimate " +
                                        kumo_test::num(worst) + " bohr)");
    const auto plain = chem::kumo::density_isosurface(mol);
    kumo_test::check(plain.vertices == mesh.vertices && plain.triangles == mesh.triangles,
                     name + ": density_isosurface gives the same mesh");

    // Potential at the vertices within the tolerance of the exact one.
    const chem::kumo::BasisSet basis{chem::kumo::FormattedCheckpoint(path)};
    const chem::kumo::detail::ElectrostaticPotential exact(
        basis, mol.density_matrix(chem::kumo::Spin::Total), mol.atoms());
    std::vector<double> expected(mesh.vertices.size());
    exact.evaluate(points, expected);
    const double difference = kumo_test::max_abs_diff(reference.potential, expected, name);
    kumo_test::check(difference <= 1e-8, name +
                                             ": potential within 1e-8 hartree/e of the exact "
                                             "one (largest difference " +
                                             kumo_test::num(difference) + ")");
    chem::kumo::SurfaceOptions exact_options;
    exact_options.potential_tolerance = 0.0;
    const auto exact_surface = chem::kumo::electrostatic_potential_surface(mol, exact_options);
    kumo_test::check(kumo_test::max_abs_diff(exact_surface.potential, expected, name) <= 1e-12,
                     name + ": tolerance 0 gives the exact potential to 1e-12");
    const auto s = chem::kumo::surface_statistics(mesh, reference.potential);
    kumo_test::check(s.area == reference.statistics.area &&
                         s.total_variance == reference.statistics.total_variance,
                     name + ": statistics are surface_statistics() of the potential");

    // Identical results with 2 and 4 threads.
    for (int threads : {2, 4}) {
        kumo_test::set_threads(threads);
        const auto again = chem::kumo::electrostatic_potential_surface(mol);
        kumo_test::check(
            again.mesh.vertices == mesh.vertices && again.mesh.triangles == mesh.triangles &&
                again.potential == reference.potential,
            name + ": identical mesh and potential with " + std::to_string(threads) + " threads");
    }
    kumo_test::set_threads(1);

    // Invalid options.
    const auto with = [&](auto change) {
        chem::kumo::SurfaceOptions o;
        change(o);
        return [o, &mol] { (void)chem::kumo::electrostatic_potential_surface(mol, o); };
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    kumo_test::check_throws<std::invalid_argument>(with([](auto& o) { o.isovalue = 0; }),
                                                   name + ": isovalue 0 is rejected");
    kumo_test::check_throws<std::invalid_argument>(with([&](auto& o) { o.isovalue = nan; }),
                                                   name + ": NaN isovalue is rejected");
    kumo_test::check_throws<std::invalid_argument>(with([](auto& o) { o.spacing = -0.1; }),
                                                   name + ": negative spacing is rejected");
    kumo_test::check_throws<std::invalid_argument>(
        with([](auto& o) { o.potential_tolerance = -1e-9; }),
        name + ": negative potential tolerance is rejected");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)chem::kumo::density_isosurface(mol, 0.99e-8); },
        name + ": an isovalue below 1e-8 is rejected");
}

} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        test_sphere();
        test_octahedron();
        test_extraction_errors();
        for (int file = 1; file < argc; ++file)
            test_molecule(argv[file]);
        if (kumo_test::failures == 0) std::cout << "molecular surfaces OK\n";
    });
}
