#pragma once
#include "chem/kumo/density.hpp"
#include <cstddef>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace chem::kumo {
// Typed property keys; the first API increment deliberately supports
// density scalars and their analytic derivatives, plus MO amplitudes.
enum class DensityComponent { Total, Alpha, Beta, Spin };
enum class PropertyKind { Density, Gradient, Hessian, Laplacian, Tau, Orbital };
enum class SpinChannel { Alpha, Beta };
struct Property {
    PropertyKind kind;
    DensityComponent component = DensityComponent::Total;
    SpinChannel channel = SpinChannel::Alpha;
    std::size_t orbital_index = 0; // one-based
    friend bool operator==(const Property&, const Property&) = default;
    friend bool operator<(const Property& a, const Property& b) noexcept {
        if (a.kind != b.kind) return a.kind < b.kind;
        if (a.component != b.component) return a.component < b.component;
        if (a.channel != b.channel) return a.channel < b.channel;
        return a.orbital_index < b.orbital_index;
    }
};
struct DensityDescriptor {
    DensityComponent component;
    constexpr Property value() const { return {PropertyKind::Density, component}; }
    constexpr Property gradient() const { return {PropertyKind::Gradient, component}; }
    constexpr Property hessian() const { return {PropertyKind::Hessian, component}; }
    constexpr Property laplacian() const { return {PropertyKind::Laplacian, component}; }
    constexpr Property tau() const { return {PropertyKind::Tau, component}; }
};
namespace density {
inline constexpr DensityDescriptor total{DensityComponent::Total};
inline constexpr DensityDescriptor alpha{DensityComponent::Alpha};
inline constexpr DensityDescriptor beta{DensityComponent::Beta};
inline constexpr DensityDescriptor spin{DensityComponent::Spin};
}
struct OrbitalDescriptor {
    SpinChannel channel;
    std::size_t index;
    Property value() const { return {PropertyKind::Orbital, DensityComponent::Total, channel, index}; }
};
namespace orbital {
inline OrbitalDescriptor alpha(std::size_t index) { return {SpinChannel::Alpha, index}; }
inline OrbitalDescriptor beta(std::size_t index) { return {SpinChannel::Beta, index}; }
}
class Request {
public:
    Request& add(Property property) { properties_.push_back(property); return *this; }
    const std::vector<Property>& properties() const noexcept { return properties_; }
private:
    std::vector<Property> properties_;
};
struct EvaluationOptions {
    int threads = 0; // 0: runtime default, 1: serial
    std::size_t block_size = 128;
};
// SoA: vector fields use x,y,z planes; symmetric Hessians use xx,xy,xz,yy,yz,zz.
struct PropertyField {
    std::size_t points = 0;
    std::vector<std::vector<double>> components;
    std::span<const double> component(std::size_t index = 0) const { return components.at(index); }
};
class Results {
public:
    const PropertyField& get(Property key) const { return fields_.at(key); }
    bool contains(Property key) const noexcept { return fields_.contains(key); }
private:
    friend class Evaluator;
    std::map<Property, PropertyField> fields_;
};
class Evaluator {
public:
    explicit Evaluator(const Molecule& molecule, EvaluationOptions options = {});
    [[nodiscard]] Results evaluate(const Grid& grid, const Request& request) const;
private:
    const Molecule* molecule_;
    EvaluationOptions options_;
};
}
