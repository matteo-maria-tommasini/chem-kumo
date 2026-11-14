#pragma once
// The density properties of the kumo command line (--property, --spin,
// --backend): one table describes every property, and the functions below
// read it.
#include "chem/kumo.hpp"
#include "cli/console.hpp"
#include "cli/options.hpp"
#include <string>
#include <utility>
#include <vector>

namespace kumo_cli {
using chem::kumo::Molecule;
using chem::kumo::PropertyKind;
using chem::kumo::Spin;

/// Density components a property is defined for (as the library defines it).
enum class Components {
    Any,            ///< total, alpha, beta and spin difference
    NotDifference,  ///< total, alpha and beta (ELF, LOL)
    TotalOnly,      ///< total density only
    DifferenceOnly, ///< spin difference only (spin polarization)
};

/// How a property is computed, as named in Cube comments and table headers.
enum class Method {
    Contraction,     ///< the density contraction (AO-P or AO-MO)
    Potential,       ///< Coulomb integrals over the AO density matrix
    OrbitalEnergies, ///< SCF orbitals weighted by their energies
};

/// One --property: the library property for each component, its names and
/// where it is defined.
struct DensityProperty {
    PropertyKind kind;
    std::string name;        ///< --property name, file suffix and column prefix
    std::string description; ///< Cube comment line, with units
    std::string summary;     ///< in "Evaluating the <summary> of the ... density"
    Components components;
    Method method;
    /// Library property of a component: (descriptor.*property)().
    chem::kumo::Property (chem::kumo::DensityDescriptor::*property)() const;
};

/// Parse --property: a comma-separated list of names, duplicates ignored,
/// user order kept; throws std::invalid_argument naming the valid ones.
std::vector<DensityProperty> select_properties(const std::string& list);

/// Throw std::invalid_argument if a property is requested for a component
/// it is not defined for.
void check_components(const std::vector<std::pair<Spin, std::string>>& components,
                      const std::vector<DensityProperty>& properties);

/// The library property of `property` for a density component.
chem::kumo::Property density_property(Spin component, const DensityProperty& property);

/// "total", "alpha", "beta" or "alpha-beta difference".
std::string component_label(Spin component);

/// Table column (and PLY property) name: "<property>_<total|alpha|beta|spin>".
std::string column_name(Spin component, const DensityProperty& property);

/// How `property` was computed, for Cube comments and table headers.
std::string method_label(const DensityProperty& property, chem::kumo::Backend used,
                         double esp_tolerance);

/// "Evaluating the total electron density", or "Evaluating the gradient norm
/// and Laplacian of the total and alpha electron densities".
std::string evaluation_summary(const std::vector<std::pair<Spin, std::string>>& components,
                               const std::vector<DensityProperty>& properties);

/// Values at every point as written to files: the gradient as its norm, the
/// other properties as they are.
std::vector<double> scalar_values(const chem::kumo::PropertyField& field, PropertyKind kind);

/// Density components and properties selected with --spin and --property, and
/// the contraction selected with --backend.
struct DensitySelection {
    std::vector<std::pair<Spin, std::string>> components; ///< with Cube filename suffix
    std::vector<DensityProperty> properties;
    chem::kumo::Backend backend = chem::kumo::Backend::AUTO;
};

/// Read --spin, --property and --backend, check them, and report the
/// contraction on the console.
DensitySelection select_density(const Options& o, const Molecule& mol, Console& console);

/// Evaluator options of a selection: backend, --esp-tolerance and the
/// console's progress bar.
chem::kumo::EvaluationOptions evaluation_options(const DensitySelection& selection,
                                                 const Options& o, Console& console);
} // namespace kumo_cli
