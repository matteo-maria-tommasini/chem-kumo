#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace chem::kumo::detail {
// AO-P is authoritative. AO-MO is eligible only after validating BOTH
// total and spin density matrices against the occupied-orbital reconstruction.
enum class ElectronicBackend { AoP, AoMo };
struct BackendDecision {
    ElectronicBackend backend = ElectronicBackend::AoP;
    double total_absolute_error = 0.0;
    double total_relative_error = 0.0;
    double spin_absolute_error = 0.0;
    double spin_relative_error = 0.0;
};
inline BackendDecision validate_orbital_backend(
    const Eigen::MatrixXd& total, const Eigen::MatrixXd& spin,
    const Eigen::MatrixXd& alpha, const Eigen::MatrixXd& beta,
    std::size_t nalpha, std::size_t nbeta, bool unrestricted,
    double relative_tolerance = 1e-8, double absolute_tolerance = 1e-10) {
    if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !std::isfinite(absolute_tolerance) || absolute_tolerance < 0.0)
        throw std::invalid_argument("invalid AO-MO validation tolerance");
    BackendDecision decision;
    const auto n = total.rows();
    const Eigen::MatrixXd& beta_coeff = unrestricted ? beta : alpha;
    if (n <= 0 || total.cols() != n || spin.rows() != n || spin.cols() != n ||
        alpha.rows() != n || (unrestricted && beta.rows() != n) ||
        nalpha > static_cast<std::size_t>(alpha.cols()) ||
        nbeta > static_cast<std::size_t>(beta_coeff.cols()) ||
        (!unrestricted && nalpha != nbeta) ||
        !total.allFinite() || !spin.allFinite() || !alpha.allFinite() ||
        (unrestricted && !beta.allFinite())) return decision;
    const Eigen::MatrixXd pa = alpha.leftCols(static_cast<Eigen::Index>(nalpha)) *
                               alpha.leftCols(static_cast<Eigen::Index>(nalpha)).transpose();
    const Eigen::MatrixXd pb = beta_coeff.leftCols(static_cast<Eigen::Index>(nbeta)) *
                               beta_coeff.leftCols(static_cast<Eigen::Index>(nbeta)).transpose();
    const Eigen::MatrixXd reconstructed_total = pa + pb;
    const Eigen::MatrixXd reconstructed_spin = pa - pb;
    const double total_error = (total - reconstructed_total).norm();
    const double spin_error = (spin - reconstructed_spin).norm();
    const double total_scale = std::max(total.norm(), reconstructed_total.norm());
    const double spin_scale = std::max(spin.norm(), reconstructed_spin.norm());
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
