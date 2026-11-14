#pragma once
#include "chem/kumo/density.hpp"

#include <cstddef>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace chem::kumo {
/// Density component: total, alpha, beta, or spin difference (alpha - beta).
enum class DensityComponent { Total, Alpha, Beta, Spin };

/// Kind of property. Number of components in the result: Gradient 3 (x, y,
/// z), Hessian 6 (xx, xy, xz, yy, yz, zz), all others 1. Orbital is a signed
/// MO amplitude; the other kinds apply to a DensityComponent.
///
/// Derived kinds (see docs/API.md for the formulas and conventions):
///  - Elf: electron localization function, in [0, 1]. Total (spin-resolved
///    form, from the alpha and beta densities), Alpha or Beta; not Spin.
///  - ReducedGradient: reduced density gradient s of the total density,
///    capped at reduced_gradient_cap. Total only.
///  - SignLambda2Density: sign(lambda2) * rho of the total density, where
///    lambda2 is the middle eigenvalue of its Hessian. Total only.
///  - Lol: localized orbital locator t / (1 + t), t = D0 / tau, in [0, 1].
///    Total (spin-resolved form, as Elf), Alpha or Beta; not Spin.
///  - InteractionRegionIndicator: IRI = |grad rho| / rho^1.1 of the total
///    density, capped at interaction_region_indicator_cap. Total only.
///  - DensityOverlapRegionsIndicator: DORI = theta / (1 + theta) of the
///    total density, in [0, 1]. Total only.
///  - HamiltonianKineticEnergy: K = tau - laplacian / 4, the kinetic energy
///    density -(1/2) sum_i psi_i lap psi_i (hartree/bohr^3). Any component.
///  - PotentialEnergyDensity: V = laplacian / 4 - 2 tau (virial theorem).
///    Total only.
///  - EnergyDensity: H = tau + V = laplacian / 4 - tau = -K. Total only.
///  - SpinPolarization: (rho_alpha - rho_beta) / (rho_alpha + rho_beta), in
///    [-1, 1]; 0 where the total density is not positive. Spin only.
///  - ShannonEntropyDensity: -rho ln rho of the total density; 0 where rho is
///    not positive. Total only.
///  - FisherInformationDensity: |grad rho|^2 / rho of the total density; 0
///    where rho is not positive. Total only.
///  - AverageLocalIonizationEnergy: ALIE sum_occ n_i |e_i| |phi_i|^2 / rho,
///    hartree, over the occupied SCF orbitals. Total only.
///  - LocalElectronAffinity: LEA -sum_vir e_i |phi_i|^2 / sum_vir |phi_i|^2,
///    hartree, over all virtual orbitals. Total only.
///  - LocalElectronAttachmentEnergy: LEAE sum_vir,e_i<0 n_i e_i |phi_i|^2 / rho,
///    hartree, over the virtual orbitals with negative energy (0 everywhere
///    if there are none). Total only.
///  - LocalElectronegativity, LocalHardness: (ALIE + LEA) / 2 and
///    (ALIE - LEA) / 2, hartree. Total only.
///    These five use the SCF orbitals and their energies, whatever the
///    backend, and need orbitals that reproduce the FCHK density matrices;
///    the last four also need virtual orbitals and are not defined for a
///    restricted open-shell calculation. Each is 0 where its denominator
///    vanishes.
///  - ElectrostaticPotential: molecular electrostatic potential in hartree/e,
///    nuclear plus electronic, from the nuclear charges and the total AO
///    density matrix, within EvaluationOptions::potential_tolerance of the
///    exact Coulomb integrals. At a grid point within
///    1e-8 bohr of a nucleus, that nucleus's own singular term is left out
///    (the potential at the nucleus). Total only.
/// Requesting a derived kind for a component it is not defined for makes
/// Evaluator::evaluate() throw std::invalid_argument.
enum class PropertyKind {
    Density,
    Gradient,
    Hessian,
    Laplacian,
    Tau,
    Orbital,
    Elf,
    ReducedGradient,
    SignLambda2Density,
    ElectrostaticPotential,
    Lol,
    InteractionRegionIndicator,
    DensityOverlapRegionsIndicator,
    HamiltonianKineticEnergy,
    PotentialEnergyDensity,
    EnergyDensity,
    SpinPolarization,
    ShannonEntropyDensity,
    FisherInformationDensity,
    AverageLocalIonizationEnergy,
    LocalElectronAffinity,
    LocalElectronAttachmentEnergy,
    LocalElectronegativity,
    LocalHardness
};

/// Upper bound of the reduced density gradient returned by the Evaluator.
/// Larger values, and points where the density is not positive (where s is
/// undefined), are reported as this value; it is far above the isovalues
/// used in NCI analysis (about 0.3 to 0.6).
inline constexpr double reduced_gradient_cap = 100.0;

/// Upper bound of the interaction region indicator, as reduced_gradient_cap:
/// larger values and points where the density is not positive are reported
/// as this value. IRI isovalues are about 1.
inline constexpr double interaction_region_indicator_cap = 100.0;

/// Edge of the cubes (bohr) that group arbitrary points for the
/// electrostatic potential in Evaluator::evaluate(points): the cubes are
/// aligned with the Cartesian axes at multiples of this edge.
inline constexpr double potential_cluster_edge = 1.6;

/// MO spin channel.
enum class SpinChannel { Alpha, Beta };

/// Key of one requested property. Build it with the density:: and orbital::
/// helpers below rather than by hand.
struct Property {
    PropertyKind kind;
    DensityComponent component = DensityComponent::Total; ///< used by density kinds
    SpinChannel channel = SpinChannel::Alpha;             ///< used by Orbital
    std::size_t orbital_index = 0;                        ///< used by Orbital; one-based
    friend bool operator==(const Property&, const Property&) = default;

    friend bool operator<(const Property& a, const Property& b) noexcept {
        if (a.kind != b.kind) return a.kind < b.kind;
        if (a.component != b.component) return a.component < b.component;
        if (a.channel != b.channel) return a.channel < b.channel;
        return a.orbital_index < b.orbital_index;
    }
};

/// Properties of one density component, e.g. density::total.gradient().
struct DensityDescriptor {
    DensityComponent component;

    constexpr Property value() const { return {PropertyKind::Density, component}; }

    constexpr Property gradient() const { return {PropertyKind::Gradient, component}; }

    constexpr Property hessian() const { return {PropertyKind::Hessian, component}; }

    constexpr Property laplacian() const { return {PropertyKind::Laplacian, component}; }

    /// Positive kinetic-energy density.
    constexpr Property tau() const { return {PropertyKind::Tau, component}; }

    /// Electron localization function (Total, Alpha or Beta).
    constexpr Property elf() const { return {PropertyKind::Elf, component}; }

    /// Reduced density gradient (Total only).
    constexpr Property reduced_gradient() const {
        return {PropertyKind::ReducedGradient, component};
    }

    /// sign(lambda2) * rho (Total only).
    constexpr Property sign_lambda2_density() const {
        return {PropertyKind::SignLambda2Density, component};
    }

    /// Molecular electrostatic potential, nuclear plus electronic (Total only).
    constexpr Property electrostatic_potential() const {
        return {PropertyKind::ElectrostaticPotential, component};
    }

    /// Localized orbital locator (Total, Alpha or Beta).
    constexpr Property lol() const { return {PropertyKind::Lol, component}; }

    /// Interaction region indicator (Total only).
    constexpr Property interaction_region_indicator() const {
        return {PropertyKind::InteractionRegionIndicator, component};
    }

    /// Density overlap regions indicator (Total only).
    constexpr Property density_overlap_regions_indicator() const {
        return {PropertyKind::DensityOverlapRegionsIndicator, component};
    }

    /// Hamiltonian kinetic energy density K (any component).
    constexpr Property hamiltonian_kinetic_energy() const {
        return {PropertyKind::HamiltonianKineticEnergy, component};
    }

    /// Potential energy density V (Total only).
    constexpr Property potential_energy_density() const {
        return {PropertyKind::PotentialEnergyDensity, component};
    }

    /// Electronic energy density H = G + V (Total only).
    constexpr Property energy_density() const { return {PropertyKind::EnergyDensity, component}; }

    /// Spin polarization (rho_alpha - rho_beta) / rho (Spin only).
    constexpr Property spin_polarization() const {
        return {PropertyKind::SpinPolarization, component};
    }

    /// Shannon entropy density -rho ln rho (Total only).
    constexpr Property shannon_entropy_density() const {
        return {PropertyKind::ShannonEntropyDensity, component};
    }

    /// Fisher information density |grad rho|^2 / rho (Total only).
    constexpr Property fisher_information_density() const {
        return {PropertyKind::FisherInformationDensity, component};
    }

    /// Average local ionization energy (Total only).
    constexpr Property average_local_ionization_energy() const {
        return {PropertyKind::AverageLocalIonizationEnergy, component};
    }

    /// Local electron affinity (Total only).
    constexpr Property local_electron_affinity() const {
        return {PropertyKind::LocalElectronAffinity, component};
    }

    /// Local electron attachment energy (Total only).
    constexpr Property local_electron_attachment_energy() const {
        return {PropertyKind::LocalElectronAttachmentEnergy, component};
    }

    /// Local electronegativity (ALIE + LEA) / 2 (Total only).
    constexpr Property local_electronegativity() const {
        return {PropertyKind::LocalElectronegativity, component};
    }

    /// Local hardness (ALIE - LEA) / 2 (Total only).
    constexpr Property local_hardness() const { return {PropertyKind::LocalHardness, component}; }
};

/// The four density components: density::total, alpha, beta, spin.
namespace density {
inline constexpr DensityDescriptor total{DensityComponent::Total};
inline constexpr DensityDescriptor alpha{DensityComponent::Alpha};
inline constexpr DensityDescriptor beta{DensityComponent::Beta};
inline constexpr DensityDescriptor spin{DensityComponent::Spin};
} // namespace density

/// One molecular orbital, e.g. orbital::alpha(5).value().
struct OrbitalDescriptor {
    SpinChannel channel;
    std::size_t index; ///< one-based

    Property value() const {
        return {PropertyKind::Orbital, DensityComponent::Total, channel, index};
    }
};

/// MO selectors by spin channel; indices are one-based.
namespace orbital {
inline OrbitalDescriptor alpha(std::size_t index) {
    return {SpinChannel::Alpha, index};
}

inline OrbitalDescriptor beta(std::size_t index) {
    return {SpinChannel::Beta, index};
}
} // namespace orbital

/// List of properties to evaluate. Order does not matter, and a property
/// requested twice is evaluated once.
class Request {
public:
    Request& add(Property property) {
        properties_.push_back(property);
        return *this;
    }

    const std::vector<Property>& properties() const noexcept { return properties_; }

private:
    std::vector<Property> properties_;
};

/// How densities are contracted. AUTO (default) uses the occupied-MO
/// contraction (AO_MO) when the orbitals were validated against the FCHK
/// density matrices at load time, and the general density-matrix contraction
/// (AO_P) otherwise. The choice is reported by Results::backend_used().
/// AO_MO was faster in the tested workloads, but this is not guaranteed.
enum class Backend { AUTO, AO_P, AO_MO };

struct EvaluationOptions {
    /// Only 0 and 1 are accepted; larger values are rejected. The value does
    /// not currently change the number of threads: in OpenMP builds, set it
    /// with OMP_NUM_THREADS.
    int threads = 0;
    /// Grid points per evaluation block, 1 to max_block_size.
    std::size_t block_size = 128;
    Backend backend = Backend::AUTO;
    /// AO screening threshold for densities and density derivatives, with
    /// both backends (see default_screening_threshold); must be finite and
    /// >= 0, and 0 gives the exact contraction.
    double screening_threshold = default_screening_threshold;
    /// Bound on the error of the electrostatic potential at each grid point,
    /// in hartree/e; must be finite and >= 0, and 0 evaluates every pair of
    /// Gaussian primitives exactly at every point. The grid is processed in
    /// clusters of 8 x 8 x 8 points: far groups of primitive pairs enter
    /// through multipole-to-local translations, far pairs through Taylor
    /// expansions about the cluster centre, the others exactly, and the sum
    /// of rigorous error bounds of the approximations is at most this value
    /// (see docs/API.md). The default 1e-8 is below the six significant
    /// digits of Cube files wherever |V| >= 2e-3.
    double potential_tolerance = 1e-8;
    /// Optional progress report (see ProgressCallback). `total` counts grid
    /// points over all passes of one evaluate() call: the derivative and
    /// AO-MO pass, if any, and one pass per AO-P density or orbital field.
    ProgressCallback progress = {};
};

/// Values of one property at all grid points, as one array per component
/// (structure of arrays): component(c)[point], points in Field order.
struct PropertyField {
    std::size_t points = 0;                      ///< grid points
    std::vector<std::vector<double>> components; ///< components x points

    /// Values of component `index`; throws std::out_of_range if it does not exist.
    std::span<const double> component(std::size_t index = 0) const { return components.at(index); }
};

/// Output of Evaluator::evaluate(), one PropertyField per requested property.
class Results {
public:
    /// Result for a requested property; throws std::out_of_range otherwise.
    const PropertyField& get(Property key) const { return fields_.at(key); }

    /// Whether the property was requested.
    bool contains(Property key) const noexcept { return fields_.contains(key); }

    /// Density-contraction backend actually used for this evaluation. For
    /// orbital-only or empty AUTO requests it is AO_P, as a neutral value:
    /// orbital amplitudes do not depend on this choice.
    [[nodiscard]] Backend backend_used() const noexcept { return backend_used_; }

private:
    friend class Evaluator;
    std::map<Property, PropertyField> fields_;
    Backend backend_used_ = Backend::AO_P;
};

/// Evaluates typed property requests (densities, their derivatives, MO
/// amplitudes) for one molecule on regular grids.
///
/// Lifetime: the Evaluator stores a pointer to `molecule`, not a copy. The
/// Molecule must outlive the Evaluator and every call to evaluate().
///
/// Usage:
///   Evaluator evaluator(molecule);
///   Request request;
///   request.add(density::total.value()).add(density::total.gradient());
///   const Results results = evaluator.evaluate(grid, request);
///   std::span<const double> rho = results.get(density::total.value()).component();
class Evaluator {
public:
    /// Throws std::invalid_argument for invalid options: block_size outside
    /// 1 to max_block_size, threads < 0 or > 1, a negative or nonfinite screening_threshold or
    /// potential_tolerance, an unknown backend, or AO_MO requested for a molecule whose orbitals
    /// are not AO-MO eligible.
    explicit Evaluator(const Molecule& molecule, EvaluationOptions options = {});

    /// Evaluate every property in `request` at all grid points, in Field
    /// order (third axis fastest). Throws std::invalid_argument for an invalid
    /// grid (see Grid), an orbital index of 0 (indices are one-based), or a
    /// derived property requested for a component it is not defined for.
    [[nodiscard]] Results evaluate(const Grid& grid, const Request& request) const;

    /// Evaluate every property in `request` at the columns of `points`
    /// (bohr), in column order. Values agree with grid evaluation at the same
    /// points to rounding; the electrostatic potential is evaluated in
    /// clusters of points within cubes of potential_cluster_edge bohr, fixed
    /// in space, within the same potential_tolerance. Results do not depend
    /// on the number of threads. Throws std::invalid_argument for a nonfinite
    /// point, and as evaluate(Grid) for the request.
    [[nodiscard]] Results
    evaluate(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
             const Request& request) const;

private:
    /// One evaluation on a grid or at points (defined in evaluator.cpp).
    class Run;

    const Molecule* molecule_;
    EvaluationOptions options_;
};
} // namespace chem::kumo
