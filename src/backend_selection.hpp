#pragma once
#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace chem::kumo::detail {
/// Backend chosen by validate_orbital_backend(): AoP contracts with the FCHK
/// atomic-orbital (AO) density matrix P, AoMo with the occupied molecular-
/// orbital (MO) coefficients.
enum class ElectronicBackend { AoP, AoMo };

/// Outcome of validate_orbital_backend(), with the Frobenius-norm errors of
/// the reconstructed total and spin density matrices.
struct BackendDecision {
    ElectronicBackend backend = ElectronicBackend::AoP;
    double total_absolute_error = 0.0;
    double total_relative_error = 0.0;
    double spin_absolute_error = 0.0;
    double spin_relative_error = 0.0;
};

/// Decide whether the occupied MOs may replace the FCHK density matrices.
/// The FCHK matrices are authoritative: AoMo is chosen only if BOTH the total
/// and the spin density matrix agree with their reconstruction from the
/// occupied MOs, P_alpha = C_alpha,occ C_alpha,occ^T (likewise for beta;
/// restricted wavefunctions use the alpha coefficients for both spins).
/// Inconsistent shapes or non-finite input give AoP with zero errors.
///
/// Each error is accepted if it is at most
/// absolute_tolerance + relative_tolerance * scale. The total scale is the
/// larger of the norms of the FCHK and reconstructed total matrices; the
/// spin scale is the largest of the two spin-matrix norms and the total
/// scale. The spin matrix P_alpha - P_beta is rebuilt as a difference
/// of two matrices of the size of the total one, so the rounding of the FCHK
/// coefficients (about nine significant digits) gives it an error that
/// scales with the total matrix, not with the spin matrix itself; relative
/// to a small spin matrix (a weakly spin-polarized state such as a
/// broken-symmetry singlet) that error would fail a check that the orbitals
/// actually pass.
inline BackendDecision
validate_orbital_backend(const Eigen::MatrixXd& total, const Eigen::MatrixXd& spin,
                         const Eigen::MatrixXd& alpha, const Eigen::MatrixXd& beta,
                         std::size_t nalpha, std::size_t nbeta, bool unrestricted,
                         double relative_tolerance = 1e-8, double absolute_tolerance = 1e-10) {
    if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !std::isfinite(absolute_tolerance) || absolute_tolerance < 0.0)
        throw std::invalid_argument("invalid AO-MO validation tolerance");
    BackendDecision decision;
    const auto n = total.rows();
    const Eigen::MatrixXd& beta_coeff = unrestricted ? beta : alpha;
    if (n <= 0 || total.cols() != n || spin.rows() != n || spin.cols() != n || alpha.rows() != n ||
        (unrestricted && beta.rows() != n) || nalpha > static_cast<std::size_t>(alpha.cols()) ||
        nbeta > static_cast<std::size_t>(beta_coeff.cols()) || (!unrestricted && nalpha != nbeta) ||
        !total.allFinite() || !spin.allFinite() || !alpha.allFinite() ||
        (unrestricted && !beta.allFinite()))
        return decision;
    const Eigen::MatrixXd pa = alpha.leftCols(static_cast<Eigen::Index>(nalpha)) *
                               alpha.leftCols(static_cast<Eigen::Index>(nalpha)).transpose();
    const Eigen::MatrixXd pb = beta_coeff.leftCols(static_cast<Eigen::Index>(nbeta)) *
                               beta_coeff.leftCols(static_cast<Eigen::Index>(nbeta)).transpose();
    const Eigen::MatrixXd reconstructed_total = pa + pb;
    const Eigen::MatrixXd reconstructed_spin = pa - pb;
    const double total_error = (total - reconstructed_total).norm();
    const double spin_error = (spin - reconstructed_spin).norm();
    const double total_scale = std::max(total.norm(), reconstructed_total.norm());
    const double spin_scale = std::max({spin.norm(), reconstructed_spin.norm(), total_scale});
    decision.total_absolute_error = total_error;
    decision.spin_absolute_error = spin_error;
    decision.total_relative_error = total_scale > 0.0 ? total_error / total_scale : total_error;
    decision.spin_relative_error = spin_scale > 0.0 ? spin_error / spin_scale : spin_error;
    if (total_error <= absolute_tolerance + relative_tolerance * total_scale &&
        spin_error <= absolute_tolerance + relative_tolerance * spin_scale)
        decision.backend = ElectronicBackend::AoMo;
    return decision;
}
} // namespace chem::kumo::detail
