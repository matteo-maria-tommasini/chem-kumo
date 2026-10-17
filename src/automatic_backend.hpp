#pragma once
#include "chem/kumo/evaluator.hpp"
#include "evaluation_plan.hpp"
#include <cstddef>
#include <stdexcept>

namespace chem::kumo::detail {
/*
 * Automatic density-contraction backend policy (chem-kumo 1.4.0).
 *
 * Prefer AO_MO when the occupied molecular orbitals reproduce the FCHK
 * spin-resolved density matrices within the explicitly validated tolerances.
 * For integer occupations, P_alpha = C_alpha,occ C_alpha,occ^T and
 * P_beta = C_beta,occ C_beta,occ^T. Molecule validates this equivalence
 * once during FCHK loading and caches the eligibility decision.
 *
 * Across the tested systems and workloads, AO_MO was faster than AO_P;
 * this is empirical evidence, NOT a claim of universal superiority. We
 * deliberately avoid an unvalidated AO-count/occupation threshold.
 *
 * AO_P remains essential for general FCHK density matrices (including
 * correlated densities), incomplete orbitals, unsupported occupations,
 * and any case where the orbital reconstruction is not sufficiently faithful.
 * Scientific eligibility takes precedence over measured performance.
 * Explicit AO_MO must fail for an ineligible representation, never silently
 * fall back. Explicit AO_P always uses the FCHK density matrix.
 *
 * Orbital-amplitude-only or empty requests do not involve a density
 * contraction; report AO_P as the neutral selection for those requests.
 */
inline Backend choose_backend(Backend requested, bool eligible,
                              std::size_t /*nao*/, std::size_t /*nalpha*/,
                              std::size_t /*nbeta*/, const EvaluationPlan& plan) {
    if (requested == Backend::AO_P) return Backend::AO_P;
    if (requested == Backend::AO_MO) {
        if (!eligible)
            throw std::invalid_argument("AO_MO requested for ineligible orbital representation");
        return Backend::AO_MO;
    }
    if (requested != Backend::AUTO)
        throw std::invalid_argument("unknown evaluation backend");
    if (!eligible) return Backend::AO_P;
    for (const auto& component : plan.components)
        if (component.ao_order() >= 0) return Backend::AO_MO;
    return Backend::AO_P;
}
} // namespace chem::kumo::detail
