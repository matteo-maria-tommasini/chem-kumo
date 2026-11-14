// Molecule construction from a Gaussian formatted checkpoint file: geometry,
// basis, AO density matrices, electron counts, MOs and orbital energies, and
// the load-time check of the orbital representation.

#include "backend_selection.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "density_matrix.hpp"
#include "mo_coefficients.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace chem::kumo {
namespace {
Eigen::MatrixXd unpack(const std::vector<double>& packed, std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / (n + 1) || packed.size() != n * (n + 1) / 2)
        throw std::runtime_error("invalid packed AO density matrix size");
    Eigen::MatrixXd p(n, n);
    std::size_t q = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            p(i, j) = packed[q];
            p(j, i) = packed[q];
            ++q;
        }
    if (!p.allFinite()) throw std::runtime_error("nonfinite AO density matrix");
    return p;
}

constexpr const char* alpha_count_label = "Number of alpha electrons";
constexpr const char* beta_count_label = "Number of beta electrons";

// Atoms from the FCHK geometry, in bohr.
std::vector<Atom> read_atoms(const FormattedCheckpoint& f) {
    const auto& zs = f.integers("Atomic numbers");
    const auto& xyz = f.reals("Current cartesian coordinates");
    if (zs.size() > std::numeric_limits<std::size_t>::max() / 3 || xyz.size() != 3 * zs.size() ||
        zs.empty())
        throw std::runtime_error("invalid FCHK geometry");
    const std::vector<double>* nuclear_charges = nullptr;
    if (f.contains("Nuclear charges")) {
        const auto& charges = f.reals("Nuclear charges");
        if (charges.size() != zs.size())
            throw std::runtime_error("invalid FCHK nuclear charge count");
        nuclear_charges = &charges;
    }
    std::vector<Atom> atoms;
    for (std::size_t i = 0; i < zs.size(); ++i) {
        if (zs[i] <= 0 || zs[i] > 118 || !std::isfinite(xyz[3 * i]) ||
            !std::isfinite(xyz[3 * i + 1]) || !std::isfinite(xyz[3 * i + 2]))
            throw std::runtime_error("invalid FCHK atom or coordinates");
        // FCHK atomic numbers identify the element. The nuclear charge comes from
        // the "Nuclear charges" record when present (with an effective core
        // potential, ECP, it holds the reduced charge); otherwise it is the
        // atomic number Z, which is right for all-electron calculations. The
        // effective charge of an ECP atom cannot be inferred from Z alone.
        const double charge =
            nuclear_charges == nullptr ? static_cast<double>(zs[i]) : (*nuclear_charges)[i];
        if (!std::isfinite(charge) || charge < 0.0)
            throw std::runtime_error("invalid FCHK nuclear charge");
        atoms.push_back({static_cast<int>(zs[i]), charge,
                         Eigen::Vector3d(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2])});
    }
    return atoms;
}

// Spin-difference AO density matrix. A closed-shell FCHK may omit it (zero
// matrix); an open-shell one, recognized from the electron counts, may not.
Eigen::MatrixXd read_spin_density(const FormattedCheckpoint& f, std::size_t n) {
    if (f.contains("Spin SCF Density")) return unpack(f.reals("Spin SCF Density"), n);
    if (f.contains(alpha_count_label) && f.contains(beta_count_label) &&
        f.integer(alpha_count_label) != f.integer(beta_count_label))
        throw std::runtime_error("open-shell FCHK lacks Spin SCF Density");
    return Eigen::MatrixXd::Zero(n, n);
}

// Numbers of alpha and beta electrons, each at most the basis size n.
std::pair<int, int> read_electron_counts(const FormattedCheckpoint& f, std::size_t n) {
    if (!f.contains(alpha_count_label) || !f.contains(beta_count_label))
        throw std::runtime_error("missing electron counts in FCHK");
    const auto na = f.integer(alpha_count_label);
    const auto nb = f.integer(beta_count_label);
    if (na < 0 || nb < 0 || na > static_cast<long long>(n) || nb > static_cast<long long>(n))
        throw std::runtime_error("invalid electron counts in FCHK");
    return {static_cast<int>(na), static_cast<int>(nb)};
}

// Orbital energies are optional; when both energies and MOs are present,
// there must be one energy per MO.
std::vector<double> read_orbital_energies(const FormattedCheckpoint& f, const std::string& label,
                                          const Eigen::MatrixXd& mos, const char* channel) {
    if (!f.contains(label)) return {};
    const auto& energies = f.reals(label);
    if (mos.cols() && !energies.empty() && energies.size() != static_cast<std::size_t>(mos.cols()))
        throw std::runtime_error(std::string(channel) + " orbital energy count mismatch");
    return energies;
}

// Decide once, at load time, whether the occupied molecular orbitals (MOs)
// reproduce the FCHK density matrices, which makes the occupied-orbital
// contraction (AO-MO) eligible. Missing or too few MOs leave the status
// Unavailable. An ineligible representation never prevents evaluation with
// the density-matrix contraction (AO-P).
ElectronicStructureInfo assess_orbital_representation(const Eigen::MatrixXd& total,
                                                      const Eigen::MatrixXd& spin,
                                                      const Eigen::MatrixXd& alpha_mos,
                                                      const Eigen::MatrixXd& beta_mos, int nalpha,
                                                      int nbeta, bool unrestricted, std::size_t n) {
    ElectronicStructureInfo info;
    info.occupied_alpha = static_cast<std::size_t>(nalpha);
    info.occupied_beta = static_cast<std::size_t>(nbeta);
    const bool has_alpha =
        alpha_mos.cols() >= nalpha && alpha_mos.rows() == static_cast<Eigen::Index>(n);
    const bool has_beta = !unrestricted || (beta_mos.cols() >= nbeta &&
                                            beta_mos.rows() == static_cast<Eigen::Index>(n));
    if (!has_alpha || !has_beta) return info;
    const auto decision = detail::validate_orbital_backend(
        total, spin, alpha_mos, beta_mos, info.occupied_alpha, info.occupied_beta, unrestricted,
        orbital_check_relative_tolerance, orbital_check_absolute_tolerance);
    info.total_absolute_error = decision.total_absolute_error;
    info.total_relative_error = decision.total_relative_error;
    info.spin_absolute_error = decision.spin_absolute_error;
    info.spin_relative_error = decision.spin_relative_error;
    info.orbital_status = decision.backend == detail::ElectronicBackend::AoMo
                              ? OrbitalRepresentationStatus::Valid
                              : OrbitalRepresentationStatus::Inconsistent;
    return info;
}
} // namespace

Molecule::Molecule(std::vector<Atom> a, BasisEvaluator b, Eigen::MatrixXd t, Eigen::MatrixXd s)
    : atoms_(std::move(a)), evaluator_(std::move(b)), total_(std::move(t)), spin_(std::move(s)) {}

// Records are read and validated in a fixed order (geometry, basis, density
// matrices, electron counts, MOs, energies), so a file with several defects
// always reports the first one in this order.
Molecule Molecule::from_fchk(const std::filesystem::path& path) {
    FormattedCheckpoint f(path);
    auto atoms = read_atoms(f);
    BasisEvaluator evaluator{BasisSet(f)};
    const auto n = evaluator.basis_size();
    auto total = unpack(f.reals("Total SCF Density"), n);
    auto spin = read_spin_density(f, n);
    Molecule result(std::move(atoms), std::move(evaluator), std::move(total), std::move(spin));
    std::tie(result.nalpha_, result.nbeta_) = read_electron_counts(f, n);

    // MO coefficients are optional here: densities only need the AO matrices.
    if (f.contains("Alpha MO coefficients"))
        result.alpha_mos_ = detail::read_mo_coefficients(f, "Alpha MO coefficients", n);
    if (f.contains("Beta MO coefficients"))
        result.beta_mos_ = detail::read_mo_coefficients(f, "Beta MO coefficients", n);
    result.unrestricted_ = result.beta_mos_.cols() != 0;
    if (result.nalpha_ != result.nbeta_ && !result.unrestricted_)
        throw std::runtime_error("open-shell FCHK lacks Beta MO coefficients");
    result.alpha_energies_ =
        read_orbital_energies(f, "Alpha Orbital Energies", result.alpha_mos_, "alpha");
    result.beta_energies_ =
        read_orbital_energies(f, "Beta Orbital Energies", result.beta_mos_, "beta");

    result.electronic_info_ = assess_orbital_representation(
        result.total_, result.spin_, result.alpha_mos_, result.beta_mos_, result.nalpha_,
        result.nbeta_, result.unrestricted_, n);
    return result;
}

Eigen::MatrixXd Molecule::density_matrix(DensityComponent component) const {
    Eigen::MatrixXd storage;
    return detail::select_density_matrix(total_, spin_, component, storage);
}
} // namespace chem::kumo
