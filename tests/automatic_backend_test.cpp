#include "automatic_backend.hpp"
#include <stdexcept>
#include <iostream>
using namespace chem::kumo;
int main() {
    try {
        const auto density_plan = detail::EvaluationPlan::from(Request{}.add(density::total.value()));
        const auto empty_plan = detail::EvaluationPlan::from(Request{});
        auto select = [&](Backend b, bool eligible, std::size_t nao,
                          std::size_t a, std::size_t c, const detail::EvaluationPlan& plan) {
            return detail::choose_backend(b, eligible, nao, a, c, plan);
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
        try { (void)select(Backend::AO_MO, false, 28, 8, 8, density_plan); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("explicit ineligible AO_MO was accepted");
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n'; return 1;
    }
}
