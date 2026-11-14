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
/// The orbitals are evaluated exactly (no atomic-orbital screening). Virtual
/// orbitals are
/// those above the occupied ones in each channel, among the MOs of the FCHK
/// (fewer than the basis functions if Gaussian removed linear dependencies).
struct OrbitalEnergySums {
    std::vector<double> occupied_density, occupied_energy, virtual_density, virtual_energy,
        attachment_energy;
};

class OrbitalEnergyFields {
public:
    /// `alpha` and `beta` are the MO coefficients (one MO per column);
    /// `virtuals` says whether the virtual sums are needed. Throws
    /// std::invalid_argument if the orbitals do not reproduce the FCHK density
    /// matrices or the orbital energies are missing; with `virtuals`, also for
    /// a restricted open-shell calculation (where, as in Multiwfn, the
    /// virtual-orbital fields are not defined) or if there are no virtual
    /// MOs. `what` names the quantity in the error messages.
    OrbitalEnergyFields(const Molecule& molecule, const BasisEvaluator& basis,
                        const Eigen::MatrixXd& alpha, const Eigen::MatrixXd& beta, bool virtuals,
                        const std::string& what)
        : basis_(&basis), virtuals_(virtuals) {
        // Molecule::from_fchk() already rejects open-shell files without beta
        // MOs, so this cannot happen today; the check keeps the restriction
        // explicit here in case that changes.
        if (virtuals && !molecule.unrestricted() &&
            molecule.occupied_count(Spin::Alpha) != molecule.occupied_count(Spin::Beta))
            throw std::invalid_argument(what +
                                        " is not defined for a restricted open-shell calculation");
        if (!molecule.electronic_structure_info().ao_mo_eligible())
            throw std::invalid_argument(what +
                                        " needs occupied MOs that reproduce the FCHK density "
                                        "matrices (an SCF calculation)");
        const int nalpha = molecule.occupied_count(Spin::Alpha);
        const int nbeta = molecule.occupied_count(Spin::Beta);
        const auto& alpha_energies = molecule.orbital_energies(Spin::Alpha);
        const auto& beta_energies = molecule.orbital_energies(Spin::Beta);
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
        return sums;
    }

private:
    const BasisEvaluator* basis_;
    bool virtuals_;
    Eigen::MatrixXd occupied_, virtual_, attachment_;
    Eigen::VectorXd occupied_energies_, virtual_energies_, attachment_energies_;
};
} // namespace chem::kumo::detail
