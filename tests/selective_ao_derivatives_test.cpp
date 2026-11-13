#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

namespace {
// Both finite and within 5e-11 * max(1, |x|, |y|).
void check(double x, double y, const std::string& label) {
    kumo_test::check_close(x, y, 5e-11, 5e-11, label);
}

void validate(const char* path) {
    const auto molecule = Molecule::from_fchk(path);
    Grid grid;
    grid.origin = Eigen::Vector3d(0.17, -0.31, 0.26);
    grid.steps << 0.16, 0.01, 0.0, -0.02, 0.21, 0.01, 0.01, -0.02, 0.19;
    grid.counts = Eigen::Vector3i(3, 4, 5);
    for (std::size_t block : {std::size_t{1}, std::size_t{7}, std::size_t{128}}) {
        Evaluator evaluator(molecule, EvaluationOptions{1, block});
        Request first, second;
        for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                               DensityComponent::Beta, DensityComponent::Spin}) {
            DensityDescriptor d{component};
            first.add(d.gradient()).add(d.tau());
            second.add(d.gradient()).add(d.tau()).add(d.hessian());
        }
        const auto low = evaluator.evaluate(grid, first);
        const auto full = evaluator.evaluate(grid, second);
        for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                               DensityComponent::Beta, DensityComponent::Spin}) {
            DensityDescriptor d{component};
            for (auto key : {d.gradient(), d.tau()}) {
                const auto& a = low.get(key);
                const auto& b = full.get(key);
                const std::string label =
                    std::string(path) + " block " + std::to_string(block) + " component " +
                    std::to_string(static_cast<int>(component)) + " kind " +
                    std::to_string(static_cast<int>(key.kind)) +
                    ": first-order request matches the request with the Hessian";
                for (std::size_t c = 0; c < a.components.size(); ++c)
                    for (std::size_t i = 0; i < a.component(c).size(); ++i)
                        check(a.component(c)[i], b.component(c)[i], label);
            }
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    return kumo_test::run([&] {
        validate(argv[1]);
        validate(argv[2]);
    });
}
