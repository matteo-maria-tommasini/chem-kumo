// Automatic backend selection (detail::choose_backend): AUTO picks AO_MO
// (density from molecular orbitals) when the molecule is eligible and the
// request needs a density, and AO_P (density matrix contraction) otherwise;
// an explicit backend is kept, and an explicit AO_MO for an ineligible
// molecule is rejected with std::invalid_argument. No arguments.
#include "automatic_backend.hpp"
#include <stdexcept>
#include <iostream>
using namespace chem::kumo;

int main() {
    try {
        const auto density_plan =
            detail::EvaluationPlan::from(Request{}.add(density::total.value()));
        const auto empty_plan = detail::EvaluationPlan::from(Request{});
        // The AO and occupied-orbital counts are kept in the calls below
        // for readability; they do not affect the choice.
        auto select = [&](Backend b, bool eligible, std::size_t /*nao*/, std::size_t /*nalpha*/,
                          std::size_t /*nbeta*/, const detail::EvaluationPlan& plan) {
            return detail::choose_backend(b, eligible, plan);
        };
        if (select(Backend::AUTO, false, 1410, 135, 135, density_plan) != Backend::AO_P ||
            select(Backend::AUTO, true, 1410, 135, 135, density_plan) != Backend::AO_MO ||
            select(Backend::AUTO, true, 28, 8, 8, density_plan) != Backend::AO_MO ||
            select(Backend::AUTO, true, 1410, 135, 135, empty_plan) != Backend::AO_P ||
            select(Backend::AUTO, true, 678, 135, 135, density_plan) != Backend::AO_MO ||
            select(Backend::AO_P, true, 1410, 135, 135, density_plan) != Backend::AO_P ||
            select(Backend::AO_MO, true, 28, 8, 8, density_plan) != Backend::AO_MO)
            throw std::runtime_error("automatic backend decision mismatch");
        bool rejected = false;
        try {
            (void)select(Backend::AO_MO, false, 28, 8, 8, density_plan);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        if (!rejected) throw std::runtime_error("explicit ineligible AO_MO was accepted");
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }
}
