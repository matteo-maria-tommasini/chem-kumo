#include "chem/kumo/field_sampler.hpp"

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

int main() {
    chem::kumo::Grid grid;
    grid.origin = Eigen::Vector3d(0.5, -0.7, 1.1);
    grid.steps << 0.4, 0.1, 0.0,
                  0.0, 0.5, 0.2,
                  0.0, 0.0, 0.6;
    grid.counts = Eigen::Vector3i(3, 4, 5);
    std::vector<double> values;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 5; ++k)
                values.push_back(2.0 + i - 2.0*j + 0.5*k);
    const chem::kumo::Field field(grid, {}, std::move(values));
    if ((field.grid().counts - grid.counts).norm() != 0 || !field.atoms().empty()) return 1;
    const chem::kumo::FieldSampler sampler(field);
    const Eigen::Vector3d point = grid.origin + grid.steps * Eigen::Vector3d(0.25, 1.5, 3.2);
    const auto inside = sampler.sample(point);
    if (!inside.inside || std::abs(inside.value - (2.0 + 0.25 - 3.0 + 1.6)) > 1e-12) return 2;
    const auto boundary = sampler.sample(grid.origin + grid.steps * Eigen::Vector3d(2.0, 3.0, 4.0));
    if (!boundary.inside || std::abs(boundary.value - (2.0 + 2.0 - 6.0 + 2.0)) > 1e-12) return 3;
    const auto outside = sampler.sample(grid.origin + grid.steps * Eigen::Vector3d(-0.01, 1.0, 1.0));
    if (outside.inside || outside.value != 0.0) return 4;
    try {
        (void)sampler.sample(Eigen::Vector3d(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0));
        return 5;
    } catch (const std::invalid_argument&) {}
    return 0;
}
