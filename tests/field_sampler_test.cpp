// Trilinear sampling on a non-orthogonal grid: exact for a linear field,
// inclusive at the far corner, zero outside, and strict on nonfinite input.
#include "chem/kumo/field_sampler.hpp"
#include "test_support.hpp"

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using kumo_test::check;

int main() {
    chem::kumo::Grid grid;
    grid.origin = Eigen::Vector3d(0.5, -0.7, 1.1);
    grid.steps << 0.4, 0.1, 0.0, 0.0, 0.5, 0.2, 0.0, 0.0, 0.6;
    grid.counts = Eigen::Vector3i(3, 4, 5);
    // f(i, j, k) = 2 + i - 2j + k/2 is linear, so interpolation is exact.
    std::vector<double> values;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 5; ++k)
                values.push_back(2.0 + i - 2.0 * j + 0.5 * k);
    const chem::kumo::Field field(grid, {}, std::move(values));
    check((field.grid().counts - grid.counts).norm() == 0 && field.atoms().empty(),
          "Field keeps the grid counts and an empty atom list");

    const chem::kumo::FieldSampler sampler(field);
    const auto inside = sampler.sample(grid.origin + grid.steps * Eigen::Vector3d(0.25, 1.5, 3.2));
    check(inside.inside && std::abs(inside.value - (2.0 + 0.25 - 3.0 + 1.6)) <= 1e-12,
          "an interior point interpolates the linear field exactly");

    const auto boundary = sampler.sample(grid.origin + grid.steps * Eigen::Vector3d(2.0, 3.0, 4.0));
    check(boundary.inside && std::abs(boundary.value - (2.0 + 2.0 - 6.0 + 2.0)) <= 1e-12,
          "the far corner is inside and returns the corner value");

    const auto outside =
        sampler.sample(grid.origin + grid.steps * Eigen::Vector3d(-0.01, 1.0, 1.0));
    check(!outside.inside && outside.value == 0.0,
          "a point outside the grid returns 0, inside=false");

    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            (void)sampler.sample(
                Eigen::Vector3d(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0));
        },
        "a NaN coordinate is rejected");
    return kumo_test::exit_code();
}
