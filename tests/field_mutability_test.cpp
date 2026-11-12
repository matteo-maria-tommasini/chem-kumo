#include "chem/kumo/field_sampler.hpp"
#include "test_support.hpp"

#include <Eigen/Core>

#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

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
    using kumo_test::check;
    using kumo_test::num;

    return kumo_test::run([] {
        chem::kumo::Field field = make_field();
        const auto& geometry = field.grid();
        const auto counts = geometry.counts;
        auto values = field.mutable_values();
        check(values.size() == 8, "mutable_values() spans the 8 samples");
        check(values.data() == field.values().data(),
              "mutable_values() and values() view the same storage");

        const chem::kumo::FieldSampler sampler(field);
        const Eigen::Vector3d point(0.25, 0.25, 0.25);
        const double initial = sampler.sample(point).value;
        check(std::abs(initial - 1.0) < 1e-14,
              "sampler interpolates 1 within 1e-14 (got " + num(initial, 17) + ")");

        values[0] = 9.0;
        check(field(0, 0, 0) == 9.0, "a write through mutable_values() is seen by operator()");
        check(field.values()[0] == 9.0, "a write through mutable_values() is seen by values()");
        const double at_origin = sampler.sample(Eigen::Vector3d::Zero()).value;
        check(std::abs(at_origin - 9.0) < 1e-14,
              "sampler sees the new corner value 9 within 1e-14 (got " + num(at_origin, 17) + ")");
        const double inside = sampler.sample(point).value;
        check(std::abs(inside - 4.375) < 1e-14,
              "sampler interpolates 4.375 after the write within 1e-14 (got " + num(inside, 17) +
                  ")");
        check((field.grid().counts - counts).norm() == 0, "writing samples keeps the grid counts");
        check(field.values().size() == 8, "writing samples keeps 8 samples");

        chem::kumo::Field copied = field;
        copied.mutable_values()[0] = -3.0;
        check(field.values()[0] == 9.0, "writing to a copy leaves the original unchanged");
        check(copied.values()[0] == -3.0, "writing to a copy changes the copy");
        check(copied.values().data() != field.values().data(), "a copy owns separate storage");

        chem::kumo::Field moved = std::move(copied);
        check(moved.values().size() == 8, "a moved-to field has 8 samples");
        check(moved.values()[0] == -3.0, "a moved-to field keeps the copied values");
        moved.mutable_values()[7] = 6.0;
        check(moved(1, 1, 1) == 6.0, "a write to the moved-to field is seen by operator()");
        check(field.values()[7] == 1.0,
              "a write to the moved-to field leaves the original unchanged");
    });
}
