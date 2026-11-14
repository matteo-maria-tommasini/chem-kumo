#include "cli/properties.hpp"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace kumo_cli {
namespace {
// Every --property, in the order of the help text and of error messages.
const std::vector<DensityProperty>& known_properties() {
    using D = chem::kumo::DensityDescriptor;
    using enum Components;
    using enum Method;
    static const std::vector<DensityProperty> known = {
        {PropertyKind::Density, "density", "Electron density (e/bohr^3)", "value", Any, Contraction,
         &D::value},
        {PropertyKind::Gradient, "gradient", "Electron density gradient norm |grad rho| (e/bohr^4)",
         "gradient norm", Any, Contraction, &D::gradient},
        {PropertyKind::Laplacian, "laplacian", "Electron density Laplacian (e/bohr^5)", "Laplacian",
         Any, Contraction, &D::laplacian},
        {PropertyKind::Tau, "tau", "Positive kinetic energy density tau (hartree/bohr^3)",
         "kinetic energy density tau", Any, Contraction, &D::tau},
        {PropertyKind::Elf, "elf", "Electron localization function ELF (dimensionless)", "ELF",
         NotDifference, Contraction, &D::elf},
        {PropertyKind::Lol, "lol", "Localized orbital locator LOL (dimensionless)", "LOL",
         NotDifference, Contraction, &D::lol},
        {PropertyKind::ReducedGradient, "rdg",
         "Reduced density gradient s (dimensionless, capped at 100)", "reduced density gradient",
         TotalOnly, Contraction, &D::reduced_gradient},
        {PropertyKind::SignLambda2Density, "sign-lambda2-rho",
         "sign(lambda2) * electron density (e/bohr^3)", "sign(lambda2) rho", TotalOnly, Contraction,
         &D::sign_lambda2_density},
        {PropertyKind::InteractionRegionIndicator, "iri",
         "Interaction region indicator |grad rho| / rho^1.1 (capped at 100)",
         "interaction region indicator", TotalOnly, Contraction, &D::interaction_region_indicator},
        {PropertyKind::DensityOverlapRegionsIndicator, "dori",
         "Density overlap regions indicator DORI (dimensionless)",
         "density overlap regions indicator", TotalOnly, Contraction,
         &D::density_overlap_regions_indicator},
        {PropertyKind::HamiltonianKineticEnergy, "kinetic-k",
         "Hamiltonian kinetic energy density K = tau - lap rho / 4 (hartree/bohr^3)",
         "kinetic energy density K", Any, Contraction, &D::hamiltonian_kinetic_energy},
        {PropertyKind::PotentialEnergyDensity, "potential-energy",
         "Potential energy density V = lap rho / 4 - 2 tau (hartree/bohr^3)",
         "potential energy density", TotalOnly, Contraction, &D::potential_energy_density},
        {PropertyKind::EnergyDensity, "energy-density",
         "Electronic energy density H = lap rho / 4 - tau (hartree/bohr^3)", "energy density",
         TotalOnly, Contraction, &D::energy_density},
        {PropertyKind::SpinPolarization, "spin-polarization",
         "Spin polarization (rho_alpha - rho_beta) / rho (dimensionless)", "spin polarization",
         DifferenceOnly, Contraction, &D::spin_polarization},
        {PropertyKind::ShannonEntropyDensity, "shannon",
         "Shannon entropy density -rho ln rho (e/bohr^3)", "Shannon entropy density", TotalOnly,
         Contraction, &D::shannon_entropy_density},
        {PropertyKind::FisherInformationDensity, "fisher",
         "Fisher information density |grad rho|^2 / rho (e/bohr^5)", "Fisher information density",
         TotalOnly, Contraction, &D::fisher_information_density},
        {PropertyKind::ElectrostaticPotential, "esp",
         "Molecular electrostatic potential, nuclear + electronic (hartree/e)",
         "electrostatic potential", TotalOnly, Potential, &D::electrostatic_potential},
        {PropertyKind::AverageLocalIonizationEnergy, "alie",
         "Average local ionization energy (hartree)", "average local ionization energy", TotalOnly,
         OrbitalEnergies, &D::average_local_ionization_energy},
        {PropertyKind::LocalElectronAffinity, "lea", "Local electron affinity (hartree)",
         "local electron affinity", TotalOnly, OrbitalEnergies, &D::local_electron_affinity},
        {PropertyKind::LocalElectronAttachmentEnergy, "leae",
         "Local electron attachment energy (hartree)", "local electron attachment energy",
         TotalOnly, OrbitalEnergies, &D::local_electron_attachment_energy},
        {PropertyKind::LocalElectronegativity, "local-electronegativity",
         "Local electronegativity (ALIE + LEA) / 2 (hartree)", "local electronegativity", TotalOnly,
         OrbitalEnergies, &D::local_electronegativity},
        {PropertyKind::LocalHardness, "local-hardness", "Local hardness (ALIE - LEA) / 2 (hartree)",
         "local hardness", TotalOnly, OrbitalEnergies, &D::local_hardness},
    };
    return known;
}

chem::kumo::DensityDescriptor descriptor(Spin component) {
    switch (component) {
    case Spin::Total: return chem::kumo::density::total;
    case Spin::Alpha: return chem::kumo::density::alpha;
    case Spin::Beta: return chem::kumo::density::beta;
    default: return chem::kumo::density::spin;
    }
}

std::string backend_label(chem::kumo::Backend backend) {
    return backend == chem::kumo::Backend::AO_MO ? "AO-MO contraction"
                                                 : "AO density-matrix contraction";
}

// Density contraction selected with --backend. AUTO uses the occupied MOs
// (AO_MO) when they reproduce the FCHK density matrices, checked when the
// FCHK is loaded, and the density matrices (AO_P) otherwise; AO_MO must not
// be forced on an FCHK that fails that check.
chem::kumo::Backend select_backend(const std::string& name, const Molecule& mol) {
    if (name.empty() || name == "auto") return chem::kumo::Backend::AUTO;
    if (name == "ao-p") return chem::kumo::Backend::AO_P;
    if (!mol.electronic_structure_info().ao_mo_eligible())
        throw std::invalid_argument("--backend ao-mo: the occupied orbitals in this FCHK do not "
                                    "reproduce its density matrices; use auto or ao-p");
    return chem::kumo::Backend::AO_MO;
}

// Which contraction --backend resolves to, and why.
std::string backend_summary(chem::kumo::Backend backend, const Molecule& mol) {
    const bool eligible = mol.electronic_structure_info().ao_mo_eligible();
    switch (backend) {
    case chem::kumo::Backend::AO_P:
        return "Density contraction: AO density matrices (requested with --backend ao-p)";
    case chem::kumo::Backend::AO_MO:
        return "Density contraction: occupied molecular orbitals (requested with --backend ao-mo)";
    default:
        return eligible ? "Density contraction: occupied molecular orbitals (AO-MO); they "
                          "reproduce the FCHK density matrices"
                        : "Density contraction: AO density matrices (AO-P); the occupied "
                          "orbitals do not reproduce the FCHK density matrices";
    }
}
} // namespace

std::vector<DensityProperty> select_properties(const std::string& list) {
    const auto& known = known_properties();
    std::vector<DensityProperty> result;
    std::size_t start = 0;
    while (true) {
        const auto comma = list.find(',', start);
        const auto token =
            list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        const auto match = std::find_if(known.begin(), known.end(),
                                        [&](const DensityProperty& p) { return p.name == token; });
        if (match == known.end()) {
            std::string names;
            for (const auto& p : known)
                names += (names.empty() ? "" : ", ") + p.name;
            throw std::invalid_argument("unknown --property '" + token + "' (use one of " + names +
                                        ")");
        }
        if (std::none_of(result.begin(), result.end(),
                         [&](const DensityProperty& p) { return p.kind == match->kind; }))
            result.push_back(*match);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return result;
}

void check_components(const std::vector<std::pair<Spin, std::string>>& components,
                      const std::vector<DensityProperty>& properties) {
    for (const auto& property : properties)
        for (const auto& component : components) {
            const Spin spin = component.first;
            if (property.components == Components::NotDifference && spin == Spin::Difference)
                throw std::invalid_argument("--property " + property.name +
                                            " is not defined for --spin difference (or all); "
                                            "use total, alpha or beta");
            if (property.components == Components::DifferenceOnly && spin != Spin::Difference)
                throw std::invalid_argument("--property " + property.name +
                                            " is defined only for --spin difference");
            if (property.components == Components::TotalOnly && spin != Spin::Total)
                throw std::invalid_argument("--property " + property.name +
                                            " is defined only for --spin total");
        }
}

chem::kumo::Property density_property(Spin component, const DensityProperty& property) {
    const auto d = descriptor(component);
    return (d.*property.property)();
}

std::string component_label(Spin component) {
    switch (component) {
    case Spin::Total: return "total";
    case Spin::Alpha: return "alpha";
    case Spin::Beta: return "beta";
    default: return "alpha-beta difference";
    }
}

std::string column_name(Spin component, const DensityProperty& property) {
    const char* names[] = {"total", "alpha", "beta", "spin"};
    const auto index = component == Spin::Total   ? 0
                       : component == Spin::Alpha ? 1
                       : component == Spin::Beta  ? 2
                                                  : 3;
    return property.name + "_" + names[index];
}

std::string method_label(const DensityProperty& property, chem::kumo::Backend used,
                         double esp_tolerance) {
    switch (property.method) {
    case Method::OrbitalEnergies: return "SCF orbitals and orbital energies";
    case Method::Potential: {
        if (esp_tolerance == 0.0) return "exact Coulomb integrals over the AO density matrix";
        std::ostringstream bound;
        bound << "error <= " << esp_tolerance << " hartree/e";
        return "Coulomb integrals over the AO density matrix, " + bound.str();
    }
    default: return backend_label(used);
    }
}

std::string evaluation_summary(const std::vector<std::pair<Spin, std::string>>& components,
                               const std::vector<DensityProperty>& properties) {
    std::vector<std::string> what, which;
    for (const auto& property : properties)
        what.push_back(property.summary);
    for (const auto& component : components)
        which.push_back(component_label(component.first));
    const auto densities = join(which) + " electron densit" + (which.size() > 1 ? "ies" : "y");
    if (what.size() == 1 && what.front() == "value") return "Evaluating the " + densities;
    return "Evaluating the " + join(what) + " of the " + densities;
}

std::vector<double> scalar_values(const chem::kumo::PropertyField& field, PropertyKind kind) {
    std::vector<double> values(field.points);
    if (kind == PropertyKind::Gradient) {
        for (std::size_t q = 0; q < field.points; ++q)
            values[q] = std::sqrt(field.component(0)[q] * field.component(0)[q] +
                                  field.component(1)[q] * field.component(1)[q] +
                                  field.component(2)[q] * field.component(2)[q]);
    } else {
        const auto source = field.component();
        std::copy(source.begin(), source.end(), values.begin());
    }
    return values;
}

DensitySelection select_density(const Options& o, const Molecule& mol, Console& console) {
    DensitySelection selection;
    auto& components = selection.components;
    const auto spin = o.spin.empty() ? "total" : o.spin;
    if (spin == "total" || spin == "all") components.emplace_back(Spin::Total, "density");
    if (spin == "alpha" || spin == "all") components.emplace_back(Spin::Alpha, "alpha");
    if (spin == "beta" || spin == "all") components.emplace_back(Spin::Beta, "beta");
    if (spin == "difference" || spin == "all") components.emplace_back(Spin::Difference, "spin");
    if (components.empty())
        throw std::invalid_argument(
            "density --spin must be total, alpha, beta, difference, or all");
    selection.properties = select_properties(o.properties.empty() ? "density" : o.properties);
    check_components(components, selection.properties);
    if (o.esp_tolerance_set &&
        std::none_of(selection.properties.begin(), selection.properties.end(),
                     [](const DensityProperty& p) {
                         return p.kind == PropertyKind::ElectrostaticPotential;
                     }))
        throw std::invalid_argument("--esp-tolerance applies only to --property esp");
    selection.backend = select_backend(o.backend, mol);
    console.say(backend_summary(selection.backend, mol));
    return selection;
}

chem::kumo::EvaluationOptions evaluation_options(const DensitySelection& selection,
                                                 const Options& o, Console& console) {
    chem::kumo::EvaluationOptions options;
    options.backend = selection.backend;
    options.potential_tolerance = o.esp_tolerance;
    options.progress = console.progress();
    return options;
}
} // namespace kumo_cli
