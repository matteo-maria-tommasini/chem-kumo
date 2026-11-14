#pragma once
#include "chem/kumo/evaluator.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

namespace chem::kumo::detail {
/// Per-component requirements of a Request: which properties are needed and,
/// from them, the highest AO derivative order to evaluate (-1: none).
///
/// Derived properties are computed pointwise from density derivatives:
/// ELF of a spin channel from its value, gradient and tau; the total ELF
/// from both channels (so it marks Alpha and Beta with `elf_input` and needs
/// no Total derivatives itself); the reduced gradient from the total value
/// and gradient; sign(lambda2) rho from the total value and Hessian. LOL
/// follows ELF (value and tau per channel, `lol_input`). IRI, the Shannon
/// and Fisher densities need the total value and gradient, DORI also the
/// Hessian, and K, V and H the Laplacian and tau of their component. The
/// spin polarization needs the Total and Spin values (`value_input`).
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
    bool lol = false;              ///< LOL of this component is requested
    bool lol_input = false;        ///< value and tau feed a LOL
    bool gradient_input = false;   ///< value and gradient feed IRI, Shannon, Fisher
    bool hessian_input = false;    ///< value, gradient and Hessian feed DORI
    bool energy_input = false;     ///< Laplacian and tau feed K, V or H
    bool value_input = false;      ///< the value feeds the spin polarization
    /// Total only; computed from the density matrix by its own pass, so it
    /// needs no AO values or derivatives.
    bool electrostatic_potential = false;

    /// Whether this component's derivatives must be evaluated on points.
    /// A total ELF alone does not need them (it uses Alpha and Beta).
    [[nodiscard]] bool derivatives() const noexcept {
        return gradient || hessian || laplacian || tau || elf_input || reduced_gradient ||
               sign_lambda2 || lol_input || gradient_input || hessian_input || energy_input ||
               value_input;
    }

    [[nodiscard]] int ao_order() const noexcept {
        if (hessian || laplacian || sign_lambda2 || hessian_input || energy_input) return 2;
        // value_input needs only values, but the AO-P point kernel starts at
        // first derivatives.
        if (gradient || tau || elf_input || reduced_gradient || lol_input || gradient_input ||
            value_input)
            return 1;
        return density ? 0 : -1;
    }
};

struct EvaluationPlan {
    std::array<ComponentPlan, 4> components{};
    bool orbitals = false;
    /// ALIE, LEA, LEAE, local electronegativity or hardness requested: they
    /// come from the SCF orbitals and energies by their own pass.
    bool orbital_energies = false;
    bool virtual_orbitals = false; ///< also the virtual-orbital sums
    int maximum_ao_order = -1;

    [[nodiscard]] bool needs_electrostatic_potential() const noexcept {
        return components[index(DensityComponent::Total)].electrostatic_potential;
    }

    [[nodiscard]] bool needs_derivative_points() const noexcept {
        return std::any_of(components.begin(), components.end(),
                           [](const ComponentPlan& c) { return c.derivatives(); });
    }

    [[nodiscard]] static constexpr std::size_t index(DensityComponent component) noexcept {
        return static_cast<std::size_t>(component);
    }

    [[nodiscard]] static const char* orbital_energy_name(PropertyKind kind) noexcept {
        switch (kind) {
        case PropertyKind::AverageLocalIonizationEnergy:
            return "the average local ionization energy";
        case PropertyKind::LocalElectronAffinity: return "the local electron affinity";
        case PropertyKind::LocalElectronAttachmentEnergy:
            return "the local electron attachment energy";
        case PropertyKind::LocalElectronegativity: return "the local electronegativity";
        default: return "the local hardness";
        }
    }

    static void require_total(const Property& property, const char* what) {
        if (property.component != DensityComponent::Total)
            throw std::invalid_argument(std::string(what) +
                                        " is defined only for the total density");
    }

    // ELF and LOL: the total is built from the alpha and beta channels, a
    // spin channel from itself; the spin difference has no such quantity.
    void require_spin_inputs(const Property& property, bool ComponentPlan::* input,
                             const char* what) {
        if (property.component == DensityComponent::Spin)
            throw std::invalid_argument(std::string(what) +
                                        " is not defined for the spin-difference density");
        if (property.component == DensityComponent::Total) {
            components[index(DensityComponent::Alpha)].*input = true;
            components[index(DensityComponent::Beta)].*input = true;
        } else {
            components[index(property.component)].*input = true;
        }
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
                plan.require_spin_inputs(property, &ComponentPlan::elf_input, "ELF");
                break;
            case PropertyKind::ReducedGradient:
                require_total(property, "the reduced density gradient");
                c.reduced_gradient = true;
                break;
            case PropertyKind::SignLambda2Density:
                require_total(property, "sign(lambda2) rho");
                c.sign_lambda2 = true;
                break;
            case PropertyKind::ElectrostaticPotential:
                require_total(property, "the electrostatic potential");
                c.electrostatic_potential = true;
                break;
            case PropertyKind::Lol:
                c.lol = true;
                plan.require_spin_inputs(property, &ComponentPlan::lol_input, "LOL");
                break;
            case PropertyKind::InteractionRegionIndicator:
                require_total(property, "the interaction region indicator");
                c.gradient_input = true;
                break;
            case PropertyKind::ShannonEntropyDensity:
                require_total(property, "the Shannon entropy density");
                c.gradient_input = true;
                break;
            case PropertyKind::FisherInformationDensity:
                require_total(property, "the Fisher information density");
                c.gradient_input = true;
                break;
            case PropertyKind::DensityOverlapRegionsIndicator:
                require_total(property, "the density overlap regions indicator");
                c.hessian_input = true;
                break;
            case PropertyKind::HamiltonianKineticEnergy: c.energy_input = true; break;
            case PropertyKind::PotentialEnergyDensity:
                require_total(property, "the potential energy density");
                c.energy_input = true;
                break;
            case PropertyKind::EnergyDensity:
                require_total(property, "the energy density");
                c.energy_input = true;
                break;
            case PropertyKind::SpinPolarization:
                if (property.component != DensityComponent::Spin)
                    throw std::invalid_argument(
                        "the spin polarization is defined only for the spin component");
                plan.components[index(DensityComponent::Total)].value_input = true;
                c.value_input = true;
                break;
            case PropertyKind::AverageLocalIonizationEnergy:
                require_total(property, "the average local ionization energy");
                plan.orbital_energies = true;
                break;
            case PropertyKind::LocalElectronAffinity:
            case PropertyKind::LocalElectronAttachmentEnergy:
            case PropertyKind::LocalElectronegativity:
            case PropertyKind::LocalHardness:
                require_total(property, orbital_energy_name(property.kind));
                plan.orbital_energies = true;
                plan.virtual_orbitals = true;
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
