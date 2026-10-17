#include "evaluation_plan.hpp"
#include <stdexcept>
using namespace chem::kumo;
using chem::kumo::detail::EvaluationPlan;
int main() {
    const auto empty = EvaluationPlan::from(Request{});
    if (empty.maximum_ao_order != -1 || empty.needs_derivative_points()) return 1;
    const auto scalar = EvaluationPlan::from(Request{}.add(density::total.value()));
    if (scalar.maximum_ao_order != 0 || scalar.needs_derivative_points()) return 2;
    const auto first = EvaluationPlan::from(Request{}.add(density::alpha.gradient()).add(density::alpha.tau()));
    if (first.maximum_ao_order != 1 || !first.needs_derivative_points() ||
        first.components[EvaluationPlan::index(DensityComponent::Alpha)].ao_order() != 1) return 3;
    const auto mixed = EvaluationPlan::from(Request{}.add(density::spin.laplacian())
        .add(density::beta.value()).add(density::spin.laplacian())
        .add(orbital::alpha(1).value()));
    if (mixed.maximum_ao_order != 2 || !mixed.orbitals ||
        mixed.components[EvaluationPlan::index(DensityComponent::Beta)].ao_order() != 0 ||
        mixed.components[EvaluationPlan::index(DensityComponent::Spin)].ao_order() != 2) return 4;
    const auto reversed = EvaluationPlan::from(Request{}.add(orbital::alpha(1).value())
        .add(density::spin.laplacian()).add(density::beta.value()));
    if (reversed.maximum_ao_order != mixed.maximum_ao_order ||
        reversed.needs_derivative_points() != mixed.needs_derivative_points()) return 5;
    try { (void)EvaluationPlan::from(Request{}.add(orbital::alpha(0).value())); }
    catch (const std::invalid_argument&) { return 0; }
    return 6;
}
