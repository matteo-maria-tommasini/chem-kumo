#include "chem/kumo/field_integrals.hpp"
#include "chem/kumo/field_sampler.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void expect_near(double actual, double expected, double tolerance = 1e-12) {
    assert(std::isfinite(actual));
    assert(std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)));
}

template <typename Exception, typename Function>
void expect_throws(Function function) {
    bool threw = false;
    try {
        function();
    } catch (const Exception&) {
        threw = true;
    }
    assert(threw);
}

chem::kumo::Field make_field(
    const Eigen::Vector3i& counts,
    const Eigen::Matrix3d& steps,
    const Eigen::Vector3d& origin,
    double (*function)(const Eigen::Vector3d&)) {
    chem::kumo::Grid grid;
    grid.counts = counts;
    grid.steps = steps;
    grid.origin = origin;

    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(counts[0]) *
                    static_cast<std::size_t>(counts[1]) *
                    static_cast<std::size_t>(counts[2]));
    for (int i = 0; i < counts[0]; ++i)
        for (int j = 0; j < counts[1]; ++j)
            for (int k = 0; k < counts[2]; ++k) {
                const Eigen::Vector3d indices(static_cast<double>(i),
                                               static_cast<double>(j),
                                               static_cast<double>(k));
                samples.push_back(function(origin + steps * indices));
            }
    return chem::kumo::Field(grid, {}, std::move(samples));
}

double affine(const Eigen::Vector3d& r) {
    return 2.0 + 3.0 * r[0] - 2.0 * r[1] + 0.5 * r[2];
}

double constant(const Eigen::Vector3d&) {
    return 2.0;
}

double square_x(const Eigen::Vector3d& r) {
    return r[0] * r[0];
}

void test_existing_conventions_and_singleton_axes() {
    chem::kumo::Grid grid;
    grid.counts = {2, 2, 2};
    grid.steps << 1.0, 0.2, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0;
    chem::kumo::Field field(grid, {}, {1, 2, 3, 4, 5, 6, 7, 8});
    expect_near(chem::kumo::grid_cell_volume(field), 1.0);
    expect_near(chem::kumo::integrate_rectangular(field), 36.0);
    expect_near(chem::kumo::integrate_trapezoidal(field), 4.5);

    grid.counts = {1, 2, 2};
    chem::kumo::Field single(grid, {}, {1, 2, 3, 4});
    expect_near(chem::kumo::integrate_rectangular(single), 10.0);
    expect_near(chem::kumo::integrate_trapezoidal(single), 2.5);
    expect_throws<std::invalid_argument>([&] {
        (void)chem::kumo::FieldSampler(single);
    });

    for (const Eigen::Vector3i& counts : {
             Eigen::Vector3i(1, 1, 1), Eigen::Vector3i(1, 3, 4),
             Eigen::Vector3i(3, 1, 4), Eigen::Vector3i(3, 4, 1)}) {
        const auto uniform = make_field(counts, grid.steps, Eigen::Vector3d::Zero(), constant);
        const double cells = static_cast<double>(counts[0] * counts[1] * counts[2]);
        const double trapezoid_extent = static_cast<double>(
            std::max(1, counts[0] - 1) * std::max(1, counts[1] - 1) *
            std::max(1, counts[2] - 1));
        expect_near(chem::kumo::integrate_rectangular(uniform), 2.0 * cells);
        expect_near(chem::kumo::integrate_trapezoidal(uniform), 2.0 * trapezoid_extent);
    }
}

void test_oblique_affine_exactness() {
    Eigen::Matrix3d cell;
    cell << 2.0, 0.4, -0.3,
            0.0, 1.5,  0.2,
            0.0, 0.0, -0.5;
    const Eigen::Vector3d origin(0.7, -0.2, 1.3);
    const Eigen::Vector3i counts(7, 6, 5);
    const Eigen::Matrix3d steps = cell *
        Eigen::Vector3d(1.0 / (counts[0] - 1),
                        1.0 / (counts[1] - 1),
                        1.0 / (counts[2] - 1)).asDiagonal();
    const auto field = make_field(counts, steps, origin, affine);

    const double volume = std::abs(cell.determinant());
    const Eigen::Vector3d center = origin + cell * Eigen::Vector3d::Constant(0.5);
    // Tensor-product trapezoidal quadrature is exact for affine functions.
    expect_near(chem::kumo::integrate_trapezoidal(field), volume * affine(center));
    expect_near(chem::kumo::grid_cell_volume(field),
                volume / static_cast<double>((counts[0] - 1) *
                                               (counts[1] - 1) *
                                               (counts[2] - 1)));
}

void test_trapezoidal_convergence() {
    double previous_error = 0.0;
    for (const int intervals : {4, 8, 16}) {
        const Eigen::Vector3i counts(intervals + 1, 2, 2);
        Eigen::Matrix3d steps = Eigen::Matrix3d::Identity();
        steps(0, 0) = 1.0 / intervals;
        const auto field = make_field(counts, steps, Eigen::Vector3d::Zero(), square_x);
        const double error = std::abs(chem::kumo::integrate_trapezoidal(field) - 1.0 / 3.0);
        // For x^2 the composite trapezoid error is h^2/6, exactly.
        expect_near(error, 1.0 / (6.0 * intervals * intervals));
        if (previous_error > 0.0)
            expect_near(previous_error / error, 4.0);
        previous_error = error;
    }
}

void test_nonfinite_samples() {
    const Eigen::Vector3i counts(3, 3, 3);
    auto field = make_field(counts, Eigen::Matrix3d::Identity(),
                            Eigen::Vector3d::Zero(), constant);
    for (const double bad : {std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()}) {
        for (const std::size_t index : {std::size_t{0}, std::size_t{13}, std::size_t{26}}) {
            field.mutable_values()[index] = bad;
            expect_throws<std::invalid_argument>([&] {
                (void)chem::kumo::integrate_rectangular(field);
            });
            expect_throws<std::invalid_argument>([&] {
                (void)chem::kumo::integrate_trapezoidal(field);
            });
            field.mutable_values()[index] = 2.0;
        }
    }
}

void test_integral_overflow() {
    chem::kumo::Grid grid;
    grid.counts = {2, 2, 2};
    grid.steps = 1.0e100 * Eigen::Matrix3d::Identity();
    chem::kumo::Field field(grid, {}, std::vector<double>(8, 1.0e100));
    expect_throws<std::overflow_error>([&] {
        (void)chem::kumo::integrate_rectangular(field);
    });
    expect_throws<std::overflow_error>([&] {
        (void)chem::kumo::integrate_trapezoidal(field);
    });
}

} // namespace

int main() {
    test_existing_conventions_and_singleton_axes();
    test_oblique_affine_exactness();
    test_trapezoidal_convergence();
    test_nonfinite_samples();
    test_integral_overflow();
}
