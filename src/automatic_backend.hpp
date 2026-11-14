#pragma once
#include "chem/kumo/evaluator.hpp"
#include "evaluation_plan.hpp"

#include <cstddef>
#include <stdexcept>

namespace chem::kumo::detail {
/**
 * Choice of the density-contraction backend for an evaluation.
 *
 * AO_P contracts atomic-orbital (AO) values with the FCHK AO density matrix
 * P; AO_MO contracts them with the coefficients C of the occupied molecular
 * orbitals (MOs). AUTO prefers AO_MO when the occupied MOs reproduce the
 * FCHK spin-resolved density matrices within the tolerances checked by
 * validate_orbital_backend() (backend_selection.hpp). For integer
 * occupations, P_alpha = C_alpha,occ C_alpha,occ^T and
 * P_beta = C_beta,occ C_beta,occ^T. Molecule checks this equivalence once
 * while loading the FCHK file and caches the result as `eligible`.
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
 * Requests with no density component (only orbital amplitudes, or nothing
 * at all) involve no density contraction; AUTO then reports AO_P as a
 * neutral value.
 */
inline Backend choose_backend(Backend requested, bool eligible, const EvaluationPlan& plan) {
    if (requested == Backend::AO_P) return Backend::AO_P;
    if (requested == Backend::AO_MO) {
        if (!eligible)
            throw std::invalid_argument("AO_MO requested for ineligible orbital representation");
        return Backend::AO_MO;
    }
    if (requested != Backend::AUTO) throw std::invalid_argument("unknown evaluation backend");
    if (!eligible) return Backend::AO_P;
    for (const auto& component : plan.components)
        if (component.ao_order() >= 0) return Backend::AO_MO;
    return Backend::AO_P;
}
} // namespace chem::kumo::detail
