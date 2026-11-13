// Grid::line and Grid::plane: their geometry and validation, evaluation on
// them against point-by-point evaluation, and against Multiwfn reference Cube
// files of water, using a line and a plane chosen to pass through Cube points.
//
// Usage: kumo-line-plane-test WATER.fchk[.gz] H2O_LOL.cube[.gz] H2O_ELF.cube[.gz]
// (the water FCHK and Cube files of tests/data/multiwfn, see
// tests/REFERENCE_CUBES.md).
#include "chem/kumo.hpp"
#include "chem/kumo/cube.hpp"
#include "test_support.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using namespace chem::kumo;
using kumo_test::check;
using kumo_test::num;

// Same tolerance as tests/multiwfn_reference_test.cpp for LOL and ELF: the
// Cube files store six significant digits.
constexpr double reference_tolerance = 2e-5;

Eigen::Vector3d point(const Grid& grid, int i, int j, int k) {
    return grid.origin + grid.steps * Eigen::Vector3d(i, j, k);
}

void check_line_geometry() {
    const Eigen::Vector3d a(0.3, -1.0, 2.0), b(1.3, 1.0, 0.0);
    const double length = (b - a).norm(); // 3
    for (const double margin : {0.0, 0.7})
        for (const double spacing : {0.01, 0.25, 1.0, 10.0}) {
            const auto label = "line, margin " + num(margin) + ", spacing " + num(spacing);
            const Grid grid = Grid::line(a, b, margin, spacing);
            const int n = grid.counts[0];
            const double expected =
                std::max(2.0, std::ceil((length + 2 * margin) / spacing * (1 - 1e-12)) + 1);
            check(n == static_cast<int>(expected) && grid.counts[1] == 1 && grid.counts[2] == 1,
                  label + ": counts (" + std::to_string(n) + ", 1, 1), expected " +
                      num(expected, 10));
            const Eigen::Vector3d direction = (b - a) / length;
            check((point(grid, 0, 0, 0) - (a - margin * direction)).norm() < 1e-13,
                  label + ": the first point is A moved back by the margin");
            check((point(grid, n - 1, 0, 0) - (b + margin * direction)).norm() < 1e-12,
                  label + ": the last point is B moved on by the margin");
            check(grid.steps.col(0).norm() <= spacing * (1 + 1e-12),
                  label + ": the spacing does not exceed the requested one");
            check(n == 2 || grid.steps.col(0).norm() * (n - 1) / (n - 2) > spacing,
                  label + ": no fewer points would do");
            const Eigen::Matrix3d axes = grid.steps.colwise().normalized();
            check((axes.transpose() * axes - Eigen::Matrix3d::Identity()).norm() < 1e-14 &&
                      std::abs(grid.steps.col(1).norm() - 1) < 1e-14 &&
                      std::abs(grid.steps.col(2).norm() - 1) < 1e-14,
                  label + ": the axes are orthogonal, the other two of unit length");
        }
    // A length that is a multiple of the spacing up to rounding keeps that
    // multiple: 8.7 / 0.3 is 28.999999999999996 or 29.000000000000004.
    const Eigen::Vector3d start(-4.445658, 0.0, -3.100976);
    check(Grid::line(start, start + Eigen::Vector3d(8.7, 0, 0), 0.0, 0.3).counts[0] == 30,
          "a line of 8.7 bohr at spacing 0.3 has 30 points");
    check(Grid::line(start, start + Eigen::Vector3d(0, 0, 6.9), 0.0, 0.3).counts[0] == 24,
          "a line of 6.9 bohr at spacing 0.3 has 24 points");
}

void check_plane_geometry() {
    const Eigen::Vector3d a(0.2, 0.1, -0.5), b(2.2, 1.1, 0.5), c(-0.4, 2.0, 1.0);
    for (const double margin : {0.0, 1.5})
        for (const double spacing : {0.05, 0.3, 50.0}) {
            const auto label = "plane, margin " + num(margin) + ", spacing " + num(spacing);
            const Grid grid = Grid::plane(a, b, c, margin, spacing);
            const Eigen::Vector3d u = grid.steps.col(0).normalized();
            const Eigen::Vector3d v = grid.steps.col(1).normalized();
            const Eigen::Vector3d normal = grid.steps.col(2);
            check(grid.counts[0] >= 2 && grid.counts[1] >= 2 && grid.counts[2] == 1,
                  label + ": counts (>= 2, >= 2, 1)");
            check((u - (b - a).normalized()).norm() < 1e-14, label + ": u runs along B - A");
            check(std::abs(u.dot(v)) < 1e-14 && (c - a).dot(v) > 0,
                  label + ": v is orthogonal to u and points towards C");
            check((normal - u.cross(v)).norm() < 1e-14, label + ": the third axis is u x v");
            for (const Eigen::Vector3d& p : {a, b, c})
                check(std::abs((p - grid.origin).dot(normal)) < 1e-13,
                      label + ": A, B and C lie in the plane");
            // In-plane extent: the three points enlarged by the margin, both
            // edges included.
            double u_low = 0, u_high = 0, v_low = 0, v_high = 0;
            for (const Eigen::Vector3d& p : {a, b, c}) {
                u_low = std::min(u_low, (p - a).dot(u));
                u_high = std::max(u_high, (p - a).dot(u));
                v_low = std::min(v_low, (p - a).dot(v));
                v_high = std::max(v_high, (p - a).dot(v));
            }
            const Eigen::Vector3d last = point(grid, grid.counts[0] - 1, grid.counts[1] - 1, 0);
            check(std::abs((grid.origin - a).dot(u) - (u_low - margin)) < 1e-13 &&
                      std::abs((grid.origin - a).dot(v) - (v_low - margin)) < 1e-13 &&
                      std::abs((last - a).dot(u) - (u_high + margin)) < 1e-12 &&
                      std::abs((last - a).dot(v) - (v_high + margin)) < 1e-12,
                  label + ": the edges are the points' extent enlarged by the margin");
            check(grid.steps.col(0).norm() <= spacing * (1 + 1e-12) &&
                      grid.steps.col(1).norm() <= spacing * (1 + 1e-12),
                  label + ": the spacings do not exceed the requested one");
        }
}

void check_rejections() {
    const Eigen::Vector3d a(0, 0, 0), b(1, 0, 0), c(0, 1, 0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const Eigen::Vector3d bad(nan, 0, 0);
    auto rejected = [](const std::string& label, auto call) {
        kumo_test::check_throws<std::invalid_argument>(call, label + " is rejected");
    };
    rejected("a line with coincident ends", [&] { (void)Grid::line(a, a, 0.0, 0.1); });
    rejected("a line with a NaN end", [&] { (void)Grid::line(a, bad, 0.0, 0.1); });
    rejected("a line with a negative margin", [&] { (void)Grid::line(a, b, -1.0, 0.1); });
    rejected("a line with an infinite margin", [&] { (void)Grid::line(a, b, inf, 0.1); });
    rejected("a line with zero spacing", [&] { (void)Grid::line(a, b, 0.0, 0.0); });
    rejected("a line with NaN spacing", [&] { (void)Grid::line(a, b, 0.0, nan); });
    rejected("a line with too many points", [&] { (void)Grid::line(a, b, 0.0, 1e-300); });
    rejected("a plane of collinear points",
             [&] { (void)Grid::plane(a, b, Eigen::Vector3d(2, 0, 0), 0.0, 0.1); });
    rejected("a plane of nearly collinear points",
             [&] { (void)Grid::plane(a, b, Eigen::Vector3d(2, 1e-10, 0), 0.0, 0.1); });
    rejected("a plane with coincident points", [&] { (void)Grid::plane(a, a, c, 0.0, 0.1); });
    rejected("a plane with a NaN point", [&] { (void)Grid::plane(a, b, bad, 0.0, 0.1); });
    rejected("a plane with a negative margin", [&] { (void)Grid::plane(a, b, c, -1.0, 0.1); });
    rejected("a plane with zero spacing", [&] { (void)Grid::plane(a, b, c, 0.0, 0.0); });
    rejected("a plane with too many points", [&] { (void)Grid::plane(a, b, c, 0.0, 1e-300); });
}

// Values on a line and a plane equal those evaluated one point at a time.
void check_against_points(const Molecule& molecule) {
    const auto& atoms = molecule.atoms();
    Request request;
    request.add(density::total.value())
        .add(density::total.laplacian())
        .add(density::total.elf())
        .add(density::total.density_overlap_regions_indicator())
        .add(density::total.electrostatic_potential());
    const Evaluator evaluator(molecule);
    const Grid grids[] = {
        Grid::line(atoms[0].position, atoms[1].position, 1.0, 0.13),
        Grid::plane(atoms[0].position, atoms[1].position, atoms[2].position, 1.0, 0.4)};
    for (const Grid& grid : grids) {
        const bool line = grid.counts[1] == 1;
        const Results results = evaluator.evaluate(grid, request);
        double worst = 0, worst_potential = 0;
        for (int i = 0; i < grid.counts[0]; ++i)
            for (int j = 0; j < grid.counts[1]; ++j) {
                Grid single;
                single.origin = point(grid, i, j, 0);
                const Results one = evaluator.evaluate(single, request);
                const auto q = static_cast<std::size_t>(i * grid.counts[1] + j);
                for (const Property& key : request.properties()) {
                    const double expected = one.get(key).component()[0];
                    const double actual = results.get(key).component()[q];
                    const double difference = std::abs(actual - expected);
                    if (key.kind == PropertyKind::ElectrostaticPotential)
                        worst_potential = std::max(worst_potential, difference);
                    else
                        worst = std::max(worst, difference / std::max(1.0, std::abs(expected)));
                }
            }
        const std::string what = line ? "line" : "plane";
        check(worst <= 1e-12,
              what + ": values equal point-by-point evaluation to 1e-12 (max " + num(worst) + ")");
        // The potential is approximated per cluster of points, within the
        // tolerance (1e-8) at each point of either evaluation.
        check(worst_potential <= 2e-8,
              what +
                  ": the electrostatic potential agrees with point-by-point evaluation "
                  "within 2e-8 (max " +
                  num(worst_potential) + ")");
    }
}

// Line and plane through Cube points: the row j = 10, k = 10 and the slice
// j = 10 (y = 0, the molecular plane of water) of the Multiwfn grid, built
// from the corner points with margin 0 and the Cube spacing.
void check_reference(const Molecule& molecule, const std::filesystem::path& cube_path, Property key,
                     const std::string& name) {
    const auto reference = read_cube(cube_path);
    const Grid& cube = reference.field.grid();
    const double spacing = cube.steps(0, 0);
    check(cube.steps.isDiagonal() && cube.steps(1, 1) == spacing && cube.steps(2, 2) == spacing,
          name + ": the reference Cube grid is cubic");
    const int nx = cube.counts[0], nz = cube.counts[2], j = 10, k = 10;

    const Grid line = Grid::line(point(cube, 0, j, k), point(cube, nx - 1, j, k), 0.0, spacing);
    const Grid plane = Grid::plane(point(cube, 0, j, 0), point(cube, nx - 1, j, 0),
                                   point(cube, 0, j, nz - 1), 0.0, spacing);
    if (!check(line.counts[0] == nx && plane.counts[0] == nx && plane.counts[1] == nz,
               name + ": the line and the plane have the points of the Cube row and slice"))
        return;

    Request request;
    request.add(key);
    const Evaluator evaluator(molecule);
    const auto compare = [&](const Grid& grid, const std::string& what) {
        const Results results = evaluator.evaluate(grid, request);
        const auto actual = results.get(key).component();
        double worst = 0, worst_position = 0;
        for (int a = 0; a < grid.counts[0]; ++a)
            for (int b = 0; b < grid.counts[1]; ++b) {
                const bool is_line = grid.counts[1] == 1;
                const int ck = is_line ? k : b;
                worst_position =
                    std::max(worst_position, (point(grid, a, b, 0) - point(cube, a, j, ck)).norm());
                const double expected = reference.field(a, j, ck);
                const double value = actual[static_cast<std::size_t>(a * grid.counts[1] + b)];
                worst = std::max(worst, std::abs(value - expected) /
                                            (reference_tolerance * std::abs(expected)));
            }
        check(worst_position < 1e-12, name + " " + what +
                                          ": the points are the Cube points (max distance " +
                                          num(worst_position) + " bohr)");
        check(worst <= 1.0, name + " " + what + ": values agree with Multiwfn within " +
                                num(reference_tolerance) + " relative (worst " + num(worst) +
                                " of the tolerance)");
        std::cout << name << " " << what << ": " << grid.counts[0] * grid.counts[1]
                  << " points, worst difference " << num(worst, 2) << " of the tolerance\n";
    };
    compare(line, "line");
    compare(plane, "plane");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 3, 3,
                                 "WATER.fchk[.gz] H2O_LOL.cube[.gz] H2O_ELF.cube[.gz]"))
        return 2;
    try {
        check_line_geometry();
        check_plane_geometry();
        check_rejections();
        const auto molecule = Molecule::from_fchk(argv[1]);
        check_against_points(molecule);
        check_reference(molecule, argv[2], density::total.lol(), "LOL");
        check_reference(molecule, argv[3], density::total.elf(), "ELF");
    } catch (const std::exception& e) {
        check(false, std::string("unexpected exception: ") + e.what());
    }
    return kumo_test::exit_code();
}
