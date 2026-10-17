#pragma once
#include "chem/kumo/evaluator.hpp"
#include <array>
#include <algorithm>

namespace chem::kumo::detail {
// Internal, immutable per-evaluation description of primitive requirements.
// This stage deliberately does not alter the validated numerical kernels.
struct ComponentPlan {
    bool density = false;
    bool gradient = false;
    bool hessian = false;
    bool laplacian = false;
    bool tau = false;
    [[nodiscard]] bool derivatives() const noexcept {
        return gradient || hessian || laplacian || tau;
    }
    [[nodiscard]] int ao_order() const noexcept {
        if (hessian || laplacian) return 2;
        if (gradient || tau) return 1;
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
                case PropertyKind::Orbital: break;
            }
        }
        for (const auto& c : plan.components)
            plan.maximum_ao_order = std::max(plan.maximum_ao_order, c.ao_order());
        return plan;
    }
};
} // namespace chem::kumo::detail
