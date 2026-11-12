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
enum class PropertyKind { Density, Gradient, Hessian, Laplacian, Tau, Orbital };

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
    std::size_t block_size = 128; ///< grid points per evaluation block; must be > 0
    Backend backend = Backend::AUTO;
    /// AO screening threshold for densities and density derivatives, with
    /// both backends (see default_screening_threshold); must be finite and
    /// >= 0, and 0 gives the exact contraction.
    double screening_threshold = default_screening_threshold;
    /// Optional progress report (see ProgressCallback). `total` counts grid
    /// points over all passes of one evaluate() call: the derivative and
    /// AO-MO pass, if any, and one pass per AO-P density or orbital field.
    ProgressCallback progress = {};
};

/// Values of one property at all grid points, as one array per component
/// (structure of arrays): component(c)[point], points in Field order.
struct PropertyField {
    std::size_t points = 0;
    std::vector<std::vector<double>> components;

    /// Values of component `index`; throws std::out_of_range if it does not exist.
    std::span<const double> component(std::size_t index = 0) const { return components.at(index); }
};

/// Output of Evaluator::evaluate(), one PropertyField per requested property.
class Results {
public:
    /// Result for a requested property; throws std::out_of_range otherwise.
    const PropertyField& get(Property key) const { return fields_.at(key); }

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
    /// Throws std::invalid_argument for invalid options: block_size == 0,
    /// threads < 0 or > 1, a negative or nonfinite screening_threshold, an
    /// unknown backend, or AO_MO requested for a molecule whose orbitals are
    /// not AO-MO eligible.
    explicit Evaluator(const Molecule& molecule, EvaluationOptions options = {});

    /// Evaluate every property in `request` at all grid points, in Field
    /// order (third axis fastest). Throws std::invalid_argument for an invalid
    /// grid (see Grid) or an orbital index of 0 (indices are one-based).
    [[nodiscard]] Results evaluate(const Grid& grid, const Request& request) const;

private:
    const Molecule* molecule_;
    EvaluationOptions options_;
};
} // namespace chem::kumo
