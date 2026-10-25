#include "chem/kumo/density.hpp"

#include <Eigen/Core>

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition) {
    if (!condition)
        throw std::runtime_error("grid geometry test failed");
}

void expect_invalid(const std::function<void()>& operation) {
    bool rejected = false;

    try {
        operation();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }

    require(rejected);
}

void test_scaled_grids() {
    for (const double scale : {1.0, 1e-6, 1e-9}) {
        chem::kumo::Grid grid;
        grid.counts = Eigen::Vector3i(3, 4, 5);
        grid.origin = Eigen::Vector3d::Zero();

        grid.steps = Eigen::Matrix3d::Identity() * scale;

        chem::kumo::Field field(
            grid, std::vector<chem::kumo::Atom>{}, std::vector<double>(static_cast<std::size_t>(grid.counts.prod()), 0.0));

        require(field.values().size() == 60);
    }
}

void test_oblique_grid() {
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(2, 3, 4);
    grid.origin = Eigen::Vector3d(1.0, 2.0, 3.0);

    grid.steps <<
        1.0, 0.2, 0.1,
        0.0, 1.0, 0.3,
        0.0, 0.0, 1.0;

    chem::kumo::Field field(
        grid, std::vector<chem::kumo::Atom>{}, std::vector<double>(static_cast<std::size_t>(grid.counts.prod()), 0.0));

    require(field.values().size() == 24);
}

void test_degenerate_grid() {
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(2, 2, 2);
    grid.origin = Eigen::Vector3d::Zero();

    grid.steps <<
        1.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 0.0;

    expect_invalid([&]() {
        chem::kumo::Field field(
            grid, std::vector<chem::kumo::Atom>{}, std::vector<double>(static_cast<std::size_t>(grid.counts.prod()), 0.0));
    });
}

void test_nonfinite_grid() {
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(2, 2, 2);
    grid.origin = Eigen::Vector3d::Zero();
    grid.steps = Eigen::Matrix3d::Identity();

    grid.origin.x() =
        std::numeric_limits<double>::quiet_NaN();

    expect_invalid([&]() {
        chem::kumo::Field field(
            grid, std::vector<chem::kumo::Atom>{}, std::vector<double>(static_cast<std::size_t>(grid.counts.prod()), 0.0));
    });
}

} // namespace

int main() {
    test_scaled_grids();
    test_oblique_grid();
    test_degenerate_grid();
    test_nonfinite_grid();
    return 0;
}
