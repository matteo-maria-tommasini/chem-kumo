// EvaluationPlan: which AO derivative order each request needs, independent
// of the order and repetition of requested properties.
#include "evaluation_plan.hpp"
#include "test_support.hpp"
#include <stdexcept>
using namespace chem::kumo;
using chem::kumo::detail::EvaluationPlan;
using kumo_test::check;

int main() {
    return kumo_test::run([&] {
        const auto empty = EvaluationPlan::from(Request{});
        check(empty.maximum_ao_order == -1 && !empty.needs_derivative_points(),
              "an empty request needs no AO evaluation");

        const auto scalar = EvaluationPlan::from(Request{}.add(density::total.value()));
        check(scalar.maximum_ao_order == 0 && !scalar.needs_derivative_points(),
              "a density value needs AO values only");

        const auto first = EvaluationPlan::from(
            Request{}.add(density::alpha.gradient()).add(density::alpha.tau()));
        check(first.maximum_ao_order == 1 && first.needs_derivative_points() &&
                  first.components[EvaluationPlan::index(DensityComponent::Alpha)].ao_order() == 1,
              "gradient and tau need first AO derivatives");

        const auto mixed = EvaluationPlan::from(Request{}
                                                    .add(density::difference.laplacian())
                                                    .add(density::beta.value())
                                                    .add(density::difference.laplacian())
                                                    .add(orbital::alpha(1).value()));
        check(
            mixed.maximum_ao_order == 2 && mixed.orbitals &&
                mixed.components[EvaluationPlan::index(DensityComponent::Beta)].ao_order() == 0 &&
                mixed.components[EvaluationPlan::index(DensityComponent::Difference)].ao_order() ==
                    2,
            "a mixed request keeps per-component orders and needs second derivatives overall");

        const auto reversed = EvaluationPlan::from(Request{}
                                                       .add(orbital::alpha(1).value())
                                                       .add(density::difference.laplacian())
                                                       .add(density::beta.value()));
        check(reversed.maximum_ao_order == mixed.maximum_ao_order &&
                  reversed.needs_derivative_points() == mixed.needs_derivative_points(),
              "the plan does not depend on request order or duplicates");

        kumo_test::check_throws<std::invalid_argument>(
            [] { (void)EvaluationPlan::from(Request{}.add(orbital::alpha(0).value())); },
            "orbital index 0 is rejected (indices are one-based)");
    });
}
