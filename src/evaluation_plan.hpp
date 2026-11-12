#pragma once
#include "chem/kumo/evaluator.hpp"
#include <array>
#include <algorithm>
#include <stdexcept>

namespace chem::kumo::detail {
/// Per-component requirements of a Request: which properties are needed and,
/// from them, the highest AO derivative order to evaluate (-1: none).
///
/// Derived properties are computed pointwise from density derivatives:
/// ELF of a spin channel from its value, gradient and tau; the total ELF
/// from both channels (so it marks Alpha and Beta with `elf_input` and needs
/// no Total derivatives itself); the reduced gradient from the total value
/// and gradient; sign(lambda2) rho from the total value and Hessian.
struct ComponentPlan {
    bool density = false;
    bool gradient = false;
    bool hessian = false;
    bool laplacian = false;
    bool tau = false;
    bool elf = false;              ///< ELF of this component is requested
    bool elf_input = false;        ///< value, gradient and tau feed an ELF
    bool reduced_gradient = false; ///< Total only
    bool sign_lambda2 = false;     ///< Total only

    /// Whether this component's derivatives must be evaluated on points.
    /// A total ELF alone does not need them (it uses Alpha and Beta).
    [[nodiscard]] bool derivatives() const noexcept {
        return gradient || hessian || laplacian || tau || elf_input || reduced_gradient ||
               sign_lambda2;
    }

    [[nodiscard]] int ao_order() const noexcept {
        if (hessian || laplacian || sign_lambda2) return 2;
        if (gradient || tau || elf_input || reduced_gradient) return 1;
        return density ? 0 : -1;
    }
};

struct EvaluationPlan {
    std::array<ComponentPlan, 4> components{};
    bool orbitals = false;
    int maximum_ao_order = -1;

    [[nodiscard]] bool needs_derivative_points() const noexcept {
        return std::any_of(components.begin(), components.end(),
                           [](const ComponentPlan& c) { return c.derivatives(); });
    }

    [[nodiscard]] static constexpr std::size_t index(DensityComponent component) noexcept {
        return static_cast<std::size_t>(component);
    }

    [[nodiscard]] static EvaluationPlan from(const Request& request) {
        EvaluationPlan plan;
        for (const Property& property : request.properties()) {
            if (property.kind == PropertyKind::Orbital) {
                if (property.orbital_index == 0)
                    throw std::invalid_argument("orbital indices are one-based");
                plan.orbitals = true;
                plan.maximum_ao_order = std::max(plan.maximum_ao_order, 0);
                continue;
            }
            auto& c = plan.components.at(index(property.component));
            switch (property.kind) {
            case PropertyKind::Density: c.density = true; break;
            case PropertyKind::Gradient: c.gradient = true; break;
            case PropertyKind::Hessian: c.hessian = true; break;
            case PropertyKind::Laplacian: c.laplacian = true; break;
            case PropertyKind::Tau: c.tau = true; break;
            case PropertyKind::Elf:
                c.elf = true;
                if (property.component == DensityComponent::Spin)
                    throw std::invalid_argument(
                        "ELF is not defined for the spin-difference density");
                if (property.component == DensityComponent::Total) {
                    plan.components[index(DensityComponent::Alpha)].elf_input = true;
                    plan.components[index(DensityComponent::Beta)].elf_input = true;
                } else {
                    c.elf_input = true;
                }
                break;
            case PropertyKind::ReducedGradient:
                if (property.component != DensityComponent::Total)
                    throw std::invalid_argument(
                        "the reduced density gradient is defined only for the total density");
                c.reduced_gradient = true;
                break;
            case PropertyKind::SignLambda2Density:
                if (property.component != DensityComponent::Total)
                    throw std::invalid_argument(
                        "sign(lambda2) rho is defined only for the total density");
                c.sign_lambda2 = true;
                break;
            case PropertyKind::Orbital: break;
            default: throw std::invalid_argument("unknown property kind");
            }
        }
        for (const auto& c : plan.components)
            plan.maximum_ao_order = std::max(plan.maximum_ao_order, c.ao_order());
        return plan;
    }
};
} // namespace chem::kumo::detail
