// Every grid-consuming entry point must apply the same geometry rules:
// positive counts, finite origin and steps, nonzero steps, and axes that are
// not (nearly) linearly dependent. The rule is scale-independent, so a valid
// grid with very small spacing is accepted everywhere.
#include "chem/kumo.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

namespace {
int failures = 0;

void fail(const std::string& message) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
}

void expect_accepted(const std::string& label, const std::function<void()>& call) {
    try {
        call();
    } catch (const std::exception& e) {
        fail(label + ": valid grid rejected (" + e.what() + ")");
    }
}

void expect_rejected(const std::string& label, const std::function<void()>& call) {
    try {
        call();
        fail(label + ": invalid grid accepted");
    } catch (const std::invalid_argument&) {
        // Expected.
    } catch (const std::exception& e) {
        fail(label + ": wrong exception type (" + e.what() + ")");
    }
}

Grid unit_grid(double spacing) {
    Grid grid;
    grid.counts = Eigen::Vector3i(2, 2, 2);
    grid.steps = spacing * Eigen::Matrix3d::Identity();
    return grid;
}

// Run the same check through every public entry point that accepts a Grid.
void check_all(const Molecule& molecule, const Grid& grid, const std::string& name, bool valid) {
    Request derivatives;
    derivatives.add(density::total.gradient());
    Request orbital_request;
    orbital_request.add(orbital::alpha(1).value());
    const Evaluator evaluator(molecule, {.threads = 1, .block_size = 4});

    const std::function<void()> calls[] = {
        [&] { (void)molecule.density(grid); },
        [&] { (void)molecule.orbital(grid, 1, Spin::Alpha); },
        [&] { (void)evaluator.evaluate(grid, derivatives); },
        [&] { (void)evaluator.evaluate(grid, orbital_request); },
    };
    const char* labels[] = {"Molecule::density", "Molecule::orbital",
                            "Evaluator (density gradient)", "Evaluator (orbital)"};
    for (std::size_t i = 0; i < std::size(calls); ++i) {
        const std::string label = name + " / " + labels[i];
        if (valid)
            expect_accepted(label, calls[i]);
        else
            expect_rejected(label, calls[i]);
    }
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " o2.fchk.gz\n";
        return 2;
    }
    const auto molecule = Molecule::from_fchk(argv[1]);

    check_all(molecule, unit_grid(0.25), "regular grid", true);
    // |det(steps)| = 1e-21: valid geometry, only the scale is small.
    check_all(molecule, unit_grid(1e-7), "small-spacing grid", true);

    Grid collinear = unit_grid(0.25);
    collinear.steps.col(2) = collinear.steps.col(0);
    check_all(molecule, collinear, "collinear axes", false);

    Grid zero_step = unit_grid(0.25);
    zero_step.steps.col(1).setZero();
    check_all(molecule, zero_step, "zero step", false);

    Grid nonfinite_origin = unit_grid(0.25);
    nonfinite_origin.origin.x() = std::numeric_limits<double>::quiet_NaN();
    check_all(molecule, nonfinite_origin, "nonfinite origin", false);

    Grid nonfinite_step = unit_grid(0.25);
    nonfinite_step.steps(0, 0) = std::numeric_limits<double>::infinity();
    check_all(molecule, nonfinite_step, "nonfinite step", false);

    Grid zero_count = unit_grid(0.25);
    zero_count.counts[1] = 0;
    check_all(molecule, zero_count, "zero count", false);

    if (failures != 0) {
        std::cerr << failures << " grid validation check(s) failed\n";
        return 1;
    }
    std::cout << "grid validation consistent across entry points\n";
    return 0;
}
