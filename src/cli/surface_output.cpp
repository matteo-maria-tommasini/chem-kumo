#include "cli/surface_output.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace kumo_cli {
using chem::kumo::Molecule;

namespace {
constexpr double bohr_in_angstrom = 0.529177210903; // CODATA 2018

// Fixed-point column, or "n/a" for an undefined value (NaN).
std::string cell(double value, int precision, int width = 14) {
    std::ostringstream out;
    if (std::isnan(value))
        out << std::setw(width) << "n/a";
    else
        out << std::fixed << std::setprecision(precision) << std::setw(width) << value;
    return out.str();
}

// Row label of a surface report.
std::string row(const std::string& name) {
    std::ostringstream out;
    out << std::left << std::setw(36) << name;
    return out.str();
}

// Area and volume table of a surface report; `areas` are (name, bohr^2).
void report_areas(const chem::kumo::TriangleMesh& mesh, double total,
                  const std::vector<std::pair<std::string, double>>& areas) {
    constexpr double a = bohr_in_angstrom;
    auto& out = std::cout;
    out << row("") << std::setw(14) << "bohr^2" << std::setw(14) << "Angstrom^2" << std::setw(10)
        << "%" << '\n';
    for (const auto& [name, value] : areas)
        out << row(name) << cell(value, 5) << cell(value * a * a, 5)
            << cell(100 * value / total, 2, 10) << '\n';
    const double volume = mesh.volume();
    out << row("Volume (bohr^3, Angstrom^3)") << cell(volume, 5) << cell(volume * a * a * a, 5)
        << "\n";
}

// Local extrema table: value in hartree(/e) and in `unit` (`factor` per
// hartree), position in Angstrom.
void report_extrema(const std::string& name, const std::vector<std::size_t>& vertices,
                    const chem::kumo::TriangleMesh& mesh, const std::vector<double>& values,
                    const std::string& hartree, const std::string& unit, double factor) {
    auto& out = std::cout;
    out << '\n'
        << name << ": " << vertices.size() << "\n"
        << "     #" << std::setw(14) << hartree << std::setw(14) << unit << std::setw(33)
        << "x, y, z (Angstrom)\n";
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        const double value = values[vertices[i]];
        const auto& x = mesh.vertices[vertices[i]];
        out << std::setw(6) << i + 1 << cell(value, 8) << cell(value * factor, 5)
            << cell(x.x() * bohr_in_angstrom, 5, 11) << cell(x.y() * bohr_in_angstrom, 5, 11)
            << cell(x.z() * bohr_in_angstrom, 5, 11) << '\n';
    }
}

// --ply: PREFIX_surface.ply with the surface and one vertex property per
// field (name, description, values); the file name follows the report.
void write_surface_ply(const Options& o, const chem::kumo::TriangleMesh& mesh,
                       const std::vector<std::pair<std::string, std::string>>& descriptions,
                       std::vector<std::pair<std::string, std::vector<double>>> fields) {
    std::ostringstream surface;
    surface << "isosurface rho = " << o.isovalue << " e/bohr^3 of the total density, grid spacing "
            << o.spacing << " bohr; x y z in bohr";
    std::vector<std::string> comments{"kumo " CHEM_KUMO_VERSION ", " + o.input.filename().string(),
                                      surface.str()};
    for (const auto& [name, description] : descriptions)
        comments.push_back(name + ": " + description);
    const auto filename = o.prefix + "_surface.ply";
    chem::kumo::write_ply(filename, mesh, fields, comments);
    std::cout << '\n' << filename << '\n';
}

// Value in a surface report, in scientific notation.
std::string scientific(double value, int width = 17) {
    std::ostringstream out;
    if (std::isnan(value))
        out << std::setw(width) << "n/a";
    else
        out << std::scientific << std::setprecision(8) << std::setw(width) << value;
    return out.str();
}
} // namespace

void report_esp_surface(const Molecule& mol, const Options& o, Console& console) {
    chem::kumo::SurfaceOptions options;
    options.isovalue = o.isovalue;
    options.spacing = o.spacing;
    options.potential_tolerance = o.esp_tolerance;
    console.begin("Computing the electrostatic potential on the density isosurface");
    const auto surface = chem::kumo::electrostatic_potential_surface(mol, options);
    console.end("Computed in");

    const auto& s = surface.statistics;
    constexpr double kcal = chem::kumo::hartree_in_kcal_per_mol;
    auto& out = std::cout;
    out << "Electrostatic potential on the total-density isosurface rho = " << o.isovalue
        << " e/bohr^3\n"
        << "Surface: " << surface.mesh.vertices.size() << " vertices, "
        << surface.mesh.triangles.size() << " triangles, grid spacing " << o.spacing << " bohr\n"
        << "Potential: Coulomb integrals over the AO density matrix, ";
    if (o.esp_tolerance > 0)
        out << "error <= " << o.esp_tolerance << " hartree/e\n\n";
    else
        out << "exact\n\n";
    out << row("") << std::setw(14) << "hartree/e" << std::setw(14) << "kcal/mol" << '\n';
    const auto potential = [&](const std::string& name, double value) {
        out << row(name) << cell(value, 8) << cell(value * kcal, 5) << '\n';
    };
    potential("Minimum", s.minimum);
    potential("Maximum", s.maximum);
    potential("Average", s.average);
    potential("Positive average", s.positive_average);
    potential("Negative average", s.negative_average);
    potential("Internal charge separation Pi", s.charge_separation);
    potential("Molecular polarity index MPI", s.polarity_index);
    out << '\n'
        << row("") << std::setw(14) << "hartree^2" << std::setw(14) << "(kcal/mol)^2" << '\n';
    const auto variance = [&](const std::string& name, double value) {
        out << row(name) << cell(value, 10) << cell(value * kcal * kcal, 5) << '\n';
    };
    variance("Total variance sigma^2_tot", s.total_variance);
    variance("Positive variance sigma^2_+", s.positive_variance);
    variance("Negative variance sigma^2_-", s.negative_variance);
    variance("nu sigma^2_tot", s.balance * s.total_variance);
    out << row("Balance of charges nu") << cell(s.balance, 8) << "\n\n";
    report_areas(surface.mesh, s.area,
                 {{"Area", s.area},
                  {"Positive area", s.positive_area},
                  {"Negative area", s.negative_area},
                  {"Nonpolar area (|V| <= 10 kcal/mol)", s.nonpolar_area}});
    report_extrema("Surface minima", s.local_minima, surface.mesh, surface.potential, "hartree/e",
                   "kcal/mol", kcal);
    report_extrema("Surface maxima", s.local_maxima, surface.mesh, surface.potential, "hartree/e",
                   "kcal/mol", kcal);
    if (o.ply)
        write_surface_ply(o, surface.mesh,
                          {{"esp", "Molecular electrostatic potential (hartree/e)"}},
                          {{"esp", surface.potential}});
}

void report_alie_surface(const Molecule& mol, const Options& o, Console& console) {
    console.begin("Computing the average local ionization energy on the density isosurface");
    const auto surface = chem::kumo::local_ionization_energy_surface(mol, o.isovalue, o.spacing);
    console.end("Computed in");

    const auto& s = surface.statistics;
    constexpr double ev = chem::kumo::hartree_in_electronvolt;
    auto& out = std::cout;
    out << "Average local ionization energy on the total-density isosurface rho = " << o.isovalue
        << " e/bohr^3\n"
        << "Surface: " << surface.mesh.vertices.size() << " vertices, "
        << surface.mesh.triangles.size() << " triangles, grid spacing " << o.spacing << " bohr\n"
        << "Ionization energy: occupied orbital energies weighted by orbital densities\n\n";
    out << row("") << std::setw(14) << "hartree" << std::setw(14) << "eV" << '\n';
    const auto energy = [&](const std::string& name, double value) {
        out << row(name) << cell(value, 8) << cell(value * ev, 5) << '\n';
    };
    energy("Minimum", s.minimum);
    energy("Maximum", s.maximum);
    energy("Average", s.average);
    energy("Average deviation Pi", s.charge_separation);
    out << '\n' << row("") << std::setw(14) << "hartree^2" << std::setw(14) << "eV^2" << '\n';
    out << row("Variance sigma^2") << cell(s.total_variance, 10)
        << cell(s.total_variance * ev * ev, 5) << "\n\n";
    report_areas(surface.mesh, s.area, {{"Area", s.area}});
    report_extrema("Surface minima", s.local_minima, surface.mesh, surface.ionization_energy,
                   "hartree", "eV", ev);
    report_extrema("Surface maxima", s.local_maxima, surface.mesh, surface.ionization_energy,
                   "hartree", "eV", ev);
    if (o.ply)
        write_surface_ply(o, surface.mesh, {{"alie", "Average local ionization energy (hartree)"}},
                          {{"alie", surface.ionization_energy}});
}

// Then (with --ply) the surface file. At most `listed` local extrema of each
// kind are printed.
void report_surface_map(const Molecule& mol, const Options& o, const DensitySelection& selection,
                        Console& console) {
    constexpr std::size_t listed = 20;
    std::vector<chem::kumo::Property> properties;
    std::vector<const DensityProperty*> mapped;              // per field
    std::vector<std::pair<std::string, std::string>> labels; // column name, description
    for (const auto& [component, suffix] : selection.components)
        for (const auto& property : selection.properties) {
            if (property.kind == PropertyKind::Gradient)
                throw std::invalid_argument(
                    "--surface-map cannot map the gradient (a vector); choose another property");
            properties.push_back(density_property(component, property));
            mapped.push_back(&property);
            // PLY property names have no hyphens.
            auto name = column_name(component, property);
            std::replace(name.begin(), name.end(), '-', '_');
            labels.emplace_back(name, property.description + ", " + component_label(component));
        }
    console.begin(evaluation_summary(selection.components, selection.properties) +
                  " on the density isosurface");
    const auto surface = chem::kumo::map_on_density_isosurface(
        mol, properties, o.isovalue, o.spacing, evaluation_options(selection, o, console));
    console.end("Computed in");

    auto& out = std::cout;
    out << "Density properties on the total-density isosurface rho = " << o.isovalue
        << " e/bohr^3\n"
        << "Surface: " << surface.mesh.vertices.size() << " vertices, "
        << surface.mesh.triangles.size() << " triangles, grid spacing " << o.spacing << " bohr\n";
    report_areas(surface.mesh, surface.statistics.front().area,
                 {{"Area", surface.statistics.front().area}});
    for (std::size_t f = 0; f < properties.size(); ++f) {
        const auto& s = surface.statistics[f];
        const auto& values = surface.values[f];
        out << "\n"
            << labels[f].first << ": " << labels[f].second << ", "
            << method_label(*mapped[f], surface.backend_used, o.esp_tolerance) << "\n\n";
        const auto line = [&](const std::string& name, double value) {
            out << row(name) << scientific(value) << '\n';
        };
        line("Minimum", s.minimum);
        line("Maximum", s.maximum);
        line("Average", s.average);
        line("Positive average", s.positive_average);
        line("Negative average", s.negative_average);
        line("Average deviation Pi", s.charge_separation);
        line("Variance sigma^2_tot", s.total_variance);
        line("Positive variance sigma^2_+", s.positive_variance);
        line("Negative variance sigma^2_-", s.negative_variance);
        out << '\n';
        report_areas(surface.mesh, s.area,
                     {{"Positive area", s.positive_area}, {"Negative area", s.negative_area}});
        const auto extrema = [&](const std::string& name, const std::vector<std::size_t>& found) {
            out << '\n'
                << name << ": " << found.size() << "\n"
                << "     #" << std::setw(17) << "value" << std::setw(33) << "x, y, z (Angstrom)\n";
            for (std::size_t i = 0; i < std::min(found.size(), listed); ++i) {
                const auto& x = surface.mesh.vertices[found[i]];
                out << std::setw(6) << i + 1 << scientific(values[found[i]])
                    << cell(x.x() * bohr_in_angstrom, 5, 11)
                    << cell(x.y() * bohr_in_angstrom, 5, 11)
                    << cell(x.z() * bohr_in_angstrom, 5, 11) << '\n';
            }
            if (found.size() > listed)
                out << "   (" << found.size() - listed << " more not listed)\n";
        };
        extrema("Surface minima", s.local_minima);
        extrema("Surface maxima", s.local_maxima);
    }
    if (o.ply) {
        std::vector<std::pair<std::string, std::vector<double>>> fields;
        for (std::size_t f = 0; f < properties.size(); ++f)
            fields.emplace_back(labels[f].first, surface.values[f]);
        write_surface_ply(o, surface.mesh, labels, std::move(fields));
    }
}
} // namespace kumo_cli
