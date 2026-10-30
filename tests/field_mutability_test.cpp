#include "chem/kumo/field_sampler.hpp"

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

void require(bool condition) {
    if (!condition)
        throw std::runtime_error("mutable field test failed");
}

chem::kumo::Field make_field() {
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(2, 2, 2);
    grid.steps = Eigen::Matrix3d::Identity();
    return chem::kumo::Field(grid, {}, std::vector<double>(8, 1.0));
}

} // namespace

int main() {
    static_assert(std::is_same_v<decltype(std::declval<chem::kumo::Field&>().mutable_values()),
                                 std::span<double>>);
    static_assert(std::is_same_v<decltype(std::declval<const chem::kumo::Field&>().values()),
                                 std::span<const double>>);

    chem::kumo::Field field = make_field();
    const auto& geometry = field.grid();
    const auto counts = geometry.counts;
    auto values = field.mutable_values();
    require(values.size() == 8);
    require(values.data() == field.values().data());

    const chem::kumo::FieldSampler sampler(field);
    const Eigen::Vector3d point(0.25, 0.25, 0.25);
    require(std::abs(sampler.sample(point).value - 1.0) < 1e-14);

    values[0] = 9.0;
    require(field(0, 0, 0) == 9.0);
    require(field.values()[0] == 9.0);
    require(std::abs(sampler.sample(Eigen::Vector3d::Zero()).value - 9.0) < 1e-14);
    require(std::abs(sampler.sample(point).value - 4.375) < 1e-14);
    require((field.grid().counts - counts).norm() == 0);
    require(field.values().size() == 8);

    chem::kumo::Field copied = field;
    copied.mutable_values()[0] = -3.0;
    require(field.values()[0] == 9.0);
    require(copied.values()[0] == -3.0);
    require(copied.values().data() != field.values().data());

    chem::kumo::Field moved = std::move(copied);
    require(moved.values().size() == 8);
    require(moved.values()[0] == -3.0);
    moved.mutable_values()[7] = 6.0;
    require(moved(1, 1, 1) == 6.0);
    require(field.values()[7] == 1.0);
    return 0;
}
