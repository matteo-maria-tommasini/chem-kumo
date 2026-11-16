#pragma once
#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/density.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace chem::kumo::detail {
/// Orbital sums behind the fields weighted by SCF orbital energies, at the
/// columns of a block of points. With phi_i the molecular orbitals (MOs), e_i
/// their energies and n_i the occupations (2 for the doubly occupied orbitals
/// of a restricted calculation, otherwise 1):
///
///   occupied_density   = sum_occ n_i |phi_i|^2              (the SCF density)
///   occupied_energy    = sum_occ n_i |e_i| |phi_i|^2
///   virtual_density    = sum_vir |phi_i|^2
///   virtual_energy     = sum_vir e_i |phi_i|^2
///   attachment_energy  = sum_(vir, e_i < 0) w e_i |phi_i|^2
///   fukui_plus         = sum_vir g_i |phi_i|^2
///   fukui_minus        = sum_occ g_i |phi_i|^2
///
/// where w = 2 in a restricted calculation (each spatial orbital counts for
/// both spins) and 1 in an unrestricted one. From these sums (Multiwfn's
/// definitions, user functions 27, -27, 28 and 29):
///
///   average local ionization energy (ALIE)
///     I(r)     = occupied_energy / occupied_density   (Sjoberg et al.)
///   local electron affinity (LEA)
///     EA_L(r)  = -virtual_energy / virtual_density    (Brinck et al.)
///   local electron attachment energy (LEAE)
///     E_att(r) = attachment_energy / occupied_density (Brinck et al. 2016)
///   local electronegativity chi_L = (I + EA_L) / 2 and
///   local hardness eta_L = (I - EA_L) / 2.
///
/// The last two sums are the orbital-weighted Fukui functions f+ and f- of
/// Pino-Rios et al. (J. Comput. Chem. 38, 481 (2017); Multiwfn's user
/// functions 95 and 96), defined for a restricted closed-shell calculation.
/// The weights g_i are Gaussians of the distance of e_i from the chemical
/// potential mu = (e_HOMO + e_LUMO) / 2, normalized over the virtual
/// orbitals (f+) or over the occupied ones (f-):
///
///   g_i = exp(-((mu - e_i) / delta)^2) / sum_j exp(-((mu - e_j) / delta)^2)
///
/// so that f+ and f- each integrate to one electron. Every orbital enters
/// the sums: unlike Multiwfn, which drops the orbitals whose unnormalized
/// Gaussian is below 1e-3, nothing is cut off.
///
/// The orbitals are evaluated without the per-block screening of atomic
/// orbitals (only the 1e-14 shell cutoff of BasisEvaluator::ao_cutoff
/// applies). Virtual
/// orbitals are those above the occupied ones in each channel, among the MOs
/// of the FCHK (fewer than the basis functions if Gaussian removed linear
/// dependencies).
struct OrbitalEnergySums {
    std::vector<double> occupied_density, occupied_energy, virtual_density, virtual_energy,
        attachment_energy, fukui_plus, fukui_minus;
};

/// Which orbital sums an OrbitalEnergyFields computes besides the occupied
/// ones (always needed).
struct OrbitalEnergyNeeds {
    /// The virtual-orbital sums (LEA, LEAE, local electronegativity and
    /// hardness).
    bool virtuals = false;
    /// The orbital-weighted Fukui sums, with Gaussian width `fukui_delta`
    /// (hartree).
    bool fukui = false;
    double fukui_delta = 0.1;
};

/// Normalized orbital-weighted Fukui weights g_i (see OrbitalEnergySums) of
/// the orbitals with the given energies. The exponents are shifted by their
/// largest value before exponentiation, which leaves the ratios unchanged
/// and keeps the largest weight away from underflow however small delta is.
inline std::vector<double> orbital_weights(const std::vector<double>& energies, double mu,
                                           double delta) {
    std::vector<double> exponents(energies.size());
    for (std::size_t i = 0; i < energies.size(); ++i) {
        const double x = (mu - energies[i]) / delta;
        exponents[i] = -x * x;
    }
    const double largest = *std::max_element(exponents.begin(), exponents.end());
    std::vector<double> weights(energies.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < energies.size(); ++i) {
        weights[i] = std::exp(exponents[i] - largest);
        sum += weights[i];
    }
    for (double& w : weights)
        w /= sum;
    return weights;
}

class OrbitalEnergyFields {
public:
    /// `alpha` and `beta` are the MO coefficients (one MO per column);
    /// `needs` says which sums are needed besides the occupied ones. Throws
    /// std::invalid_argument if the orbitals do not reproduce the FCHK density
    /// matrices or the orbital energies are missing; with `needs.virtuals`,
    /// also for a restricted open-shell calculation (where, as in Multiwfn,
    /// the virtual-orbital fields are not defined) or if there are no virtual
    /// MOs; with `needs.fukui`, also unless the calculation is restricted and
    /// closed-shell with at least one virtual MO, or if `needs.fukui_delta` is
    /// not finite and positive. `what` names the quantity in the error
    /// messages.
    OrbitalEnergyFields(const Molecule& molecule, const BasisEvaluator& basis,
                        const Eigen::MatrixXd& alpha, const Eigen::MatrixXd& beta,
                        const OrbitalEnergyNeeds& needs, const std::string& what)
        : basis_(&basis), virtuals_(needs.virtuals), fukui_(needs.fukui) {
        const bool virtuals = needs.virtuals;
        // Molecule::from_fchk() already rejects open-shell files without beta
        // MOs, so this cannot happen today; the check keeps the restriction
        // explicit here in case that changes.
        if (virtuals && !molecule.unrestricted() &&
            molecule.occupied_count(SpinChannel::Alpha) !=
                molecule.occupied_count(SpinChannel::Beta))
            throw std::invalid_argument(what +
                                        " is not defined for a restricted open-shell calculation");
        if (needs.fukui &&
            (molecule.unrestricted() || molecule.occupied_count(SpinChannel::Alpha) !=
                                            molecule.occupied_count(SpinChannel::Beta)))
            throw std::invalid_argument(
                what + " is defined only for a restricted closed-shell calculation");
        if (needs.fukui && (!std::isfinite(needs.fukui_delta) || needs.fukui_delta <= 0.0))
            throw std::invalid_argument(what + ": the Fukui delta must be finite and positive");
        if (!molecule.electronic_structure_info().ao_mo_eligible())
            throw std::invalid_argument(what +
                                        " needs occupied MOs that reproduce the FCHK density "
                                        "matrices (an SCF calculation)");
        const int nalpha = static_cast<int>(molecule.occupied_count(SpinChannel::Alpha));
        const int nbeta = static_cast<int>(molecule.occupied_count(SpinChannel::Beta));
        const auto& alpha_energies = molecule.orbital_energies(SpinChannel::Alpha);
        const auto& beta_energies = molecule.orbital_energies(SpinChannel::Beta);
        const bool unrestricted = molecule.unrestricted();
        if (alpha_energies.size() < static_cast<std::size_t>(alpha.cols()) ||
            (unrestricted && beta_energies.size() < static_cast<std::size_t>(beta.cols())) ||
            alpha_energies.size() < static_cast<std::size_t>(nalpha) ||
            beta_energies.size() < static_cast<std::size_t>(nbeta))
            throw std::invalid_argument(what + " needs the FCHK orbital energies");

        // Append one MO column and its energy. The callers scale the column by
        // sqrt(n_i) (or sqrt(w)), so that its squared amplitude carries the
        // occupation weight.
        const auto add = [](Eigen::MatrixXd& columns, Eigen::VectorXd& energies, Eigen::Index& n,
                            const Eigen::VectorXd& column, double energy) {
            columns.col(n) = column;
            energies[n] = energy;
            ++n;
        };
        const Eigen::Index rows = alpha.rows();
        occupied_.resize(rows, nalpha + nbeta);
        occupied_energies_.resize(nalpha + nbeta);
        virtual_.resize(rows, alpha.cols() + (unrestricted ? beta.cols() : 0));
        virtual_energies_.resize(virtual_.cols());
        attachment_.resize(rows, virtual_.cols());
        attachment_energies_.resize(virtual_.cols());
        Eigen::Index o = 0, v = 0, a = 0;
        const auto channel = [&](const Eigen::MatrixXd& mos, const std::vector<double>& energies,
                                 int occupied, int doubly, double weight) {
            for (Eigen::Index i = 0; i < mos.cols(); ++i) {
                const double e = energies[static_cast<std::size_t>(i)];
                if (i < occupied) {
                    const double n = i < doubly ? 2.0 : 1.0;
                    add(occupied_, occupied_energies_, o, std::sqrt(n) * mos.col(i), std::abs(e));
                } else {
                    add(virtual_, virtual_energies_, v, mos.col(i), e);
                    if (e < 0.0)
                        add(attachment_, attachment_energies_, a, std::sqrt(weight) * mos.col(i),
                            e);
                }
            }
        };
        if (unrestricted) {
            channel(alpha, alpha_energies, nalpha, 0, 1.0);
            channel(beta, beta_energies, nbeta, 0, 1.0);
        } else {
            // Restricted: the first min(nalpha, nbeta) orbitals are doubly
            // occupied, the rest up to max(nalpha, nbeta) singly occupied.
            channel(alpha, alpha_energies, std::max(nalpha, nbeta), std::min(nalpha, nbeta), 2.0);
        }
        occupied_.conservativeResize(Eigen::NoChange, o);
        occupied_energies_.conservativeResize(o);
        virtual_.conservativeResize(Eigen::NoChange, v);
        virtual_energies_.conservativeResize(v);
        attachment_.conservativeResize(Eigen::NoChange, a);
        attachment_energies_.conservativeResize(a);
        if (virtuals && v == 0)
            throw std::invalid_argument(what + " needs virtual orbitals, and the FCHK has none");
        if (needs.fukui) fukui_columns(alpha, alpha_energies, nalpha, needs.fukui_delta, what);
    }

    /// The sums at the columns of `points` (one block).
    [[nodiscard]] OrbitalEnergySums
    evaluate(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points) const {
        const Eigen::Index count = points.cols();
        Eigen::MatrixXd phi(static_cast<Eigen::Index>(basis_->basis_size()), count);
        basis_->evaluate_batch(points, phi);
        OrbitalEnergySums sums;
        const auto weighted = [&](const Eigen::MatrixXd& columns, const Eigen::VectorXd& energies,
                                  std::vector<double>& density, std::vector<double>* energy) {
            const Eigen::MatrixXd squared = (columns.transpose() * phi).array().square().matrix();
            density.resize(static_cast<std::size_t>(count));
            if (energy) energy->resize(static_cast<std::size_t>(count));
            for (Eigen::Index t = 0; t < count; ++t) {
                density[static_cast<std::size_t>(t)] = squared.col(t).sum();
                if (energy) (*energy)[static_cast<std::size_t>(t)] = energies.dot(squared.col(t));
            }
        };
        weighted(occupied_, occupied_energies_, sums.occupied_density, &sums.occupied_energy);
        if (virtuals_) {
            weighted(virtual_, virtual_energies_, sums.virtual_density, &sums.virtual_energy);
            std::vector<double> unused;
            weighted(attachment_, attachment_energies_, unused, &sums.attachment_energy);
        }
        if (fukui_) {
            weighted(fukui_plus_, Eigen::VectorXd(), sums.fukui_plus, nullptr);
            weighted(fukui_minus_, Eigen::VectorXd(), sums.fukui_minus, nullptr);
        }
        return sums;
    }

private:
    // The columns of the orbital-weighted Fukui sums of a restricted
    // closed-shell calculation (one set of spatial orbitals, the first
    // `occupied` of them doubly occupied): each MO scaled by sqrt(g_i), so
    // that its squared amplitude carries the weight. Orbitals whose weight
    // underflows to zero are left out; they would add exactly zero.
    void fukui_columns(const Eigen::MatrixXd& mos, const std::vector<double>& energies,
                       int occupied, double delta, const std::string& what) {
        const auto count = static_cast<std::size_t>(mos.cols());
        const auto homo = static_cast<std::size_t>(occupied);
        if (homo == 0 || homo >= count)
            throw std::invalid_argument(what +
                                        " needs occupied and virtual orbitals, and the FCHK lacks "
                                        "one kind");
        const double mu = 0.5 * (energies[homo - 1] + energies[homo]);
        const auto columns = [&](std::size_t first, std::size_t last, Eigen::MatrixXd& out) {
            const std::vector<double> range(energies.begin() + static_cast<std::ptrdiff_t>(first),
                                            energies.begin() + static_cast<std::ptrdiff_t>(last));
            const auto weights = orbital_weights(range, mu, delta);
            out.resize(mos.rows(), static_cast<Eigen::Index>(range.size()));
            Eigen::Index n = 0;
            for (std::size_t i = 0; i < range.size(); ++i)
                if (weights[i] > 0.0)
                    out.col(n++) =
                        std::sqrt(weights[i]) * mos.col(static_cast<Eigen::Index>(first + i));
            out.conservativeResize(Eigen::NoChange, n);
        };
        columns(homo, count, fukui_plus_);
        columns(0, homo, fukui_minus_);
    }

    const BasisEvaluator* basis_;
    bool virtuals_, fukui_;
    Eigen::MatrixXd occupied_, virtual_, attachment_, fukui_plus_, fukui_minus_;
    Eigen::VectorXd occupied_energies_, virtual_energies_, attachment_energies_;
};
} // namespace chem::kumo::detail
