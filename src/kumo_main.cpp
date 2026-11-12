#include "chem/kumo.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using chem::kumo::Molecule;
using chem::kumo::Spin;
using chem::kumo::PropertyKind;
namespace fs = std::filesystem;

namespace {
// CLI options only: numerical evaluation stays in libkumo, not in main().
struct Options {
    fs::path input;
    std::string selector;
    std::string spin;
    std::string properties;
    std::string backend;
    std::string prefix;
    double spacing = 0.20;
    double margin = 5.0;
    bool mo_mode = false;
    bool list = false;
    bool output_set = false;
};

void usage() {
    std::cout << "kumo " CHEM_KUMO_VERSION
                 " - ground-state electron densities and molecular orbitals\n\n"
              << "Usage: kumo FILE.fchk[.gz] [options]\n\n"
              << "  (default)                  Write total electron density\n"
              << "  --spin total|alpha|beta|difference|all   Density component(s)\n"
              << "  --property density,gradient,laplacian,tau\n"
              << "                             Density properties to write (default density;\n"
              << "                             gradient is the norm |grad rho|)\n"
              << "  --backend auto|ao-p|ao-mo  Density contraction (default auto: occupied MOs\n"
              << "                             when they reproduce the FCHK density matrices)\n"
              << "  --mo N[,M:K,...]           Write 1-based molecular orbitals\n"
              << "  --mo homo|lumo|frontier    Frontier orbital selection\n"
              << "  --mo homo-2:lumo+2         Frontier-relative range\n"
              << "  --spin alpha|beta|both     MO spin channel (default alpha)\n"
              << "  --list-mo                  List available orbitals and energies\n"
              << "  --spacing BOHR             Grid step (default 0.20 bohr)\n"
              << "  --margin BOHR              Bounding-box margin (default 5.0 bohr)\n"
              << "  --output PREFIX            Output filename prefix\n"
              << "  --help                     Show help\n"
              << "  --version                  Show version\n";
}

std::string lower(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

double real(const std::string& s) {
    std::size_t n = 0;
    const auto v = std::stod(s, &n);
    if (n != s.size() || !std::isfinite(v))
        throw std::invalid_argument("invalid real number: " + s);
    return v;
}

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument("missing argument after " + a);
            return argv[++i];
        };
        if (a == "--mo") {
            if (o.mo_mode) throw std::invalid_argument("--mo specified twice");
            o.selector = lower(value());
            o.mo_mode = true;
        } else if (a == "--spin") {
            if (!o.spin.empty()) throw std::invalid_argument("--spin specified twice");
            o.spin = lower(value());
        } else if (a == "--property") {
            if (!o.properties.empty()) throw std::invalid_argument("--property specified twice");
            o.properties = lower(value());
            if (o.properties.empty()) throw std::invalid_argument("empty --property list");
        } else if (a == "--backend") {
            if (!o.backend.empty()) throw std::invalid_argument("--backend specified twice");
            o.backend = lower(value());
            if (o.backend != "auto" && o.backend != "ao-p" && o.backend != "ao-mo")
                throw std::invalid_argument("--backend must be auto, ao-p or ao-mo");
        } else if (a == "--spacing")
            o.spacing = real(value());
        else if (a == "--margin")
            o.margin = real(value());
        else if (a == "--output") {
            o.prefix = value();
            o.output_set = true;
        } else if (a == "--list-mo")
            o.list = true;
        else if (!a.empty() && a[0] == '-')
            throw std::invalid_argument("unknown option: " + a);
        else if (o.input.empty())
            o.input = a;
        else
            throw std::invalid_argument("unexpected positional argument: " + a);
    }
    if (o.input.empty()) throw std::invalid_argument("missing input FCHK file");
    if (o.mo_mode && o.list) throw std::invalid_argument("choose --mo or --list-mo");
    if (!o.properties.empty() && (o.mo_mode || o.list))
        throw std::invalid_argument(
            "--property applies only to densities, not to --mo or --list-mo");
    if (!o.backend.empty() && (o.mo_mode || o.list))
        throw std::invalid_argument(
            "--backend applies only to densities, not to --mo or --list-mo");
    if (!(o.spacing > 0) || o.margin < 0)
        throw std::invalid_argument("invalid grid spacing or margin");
    if (!o.output_set) o.prefix = o.input.stem().string();
    if (o.prefix.empty()) throw std::invalid_argument("empty output prefix");
    return o;
}

std::size_t number(const std::string& s) {
    if (s.empty() ||
        !std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c); }))
        throw std::invalid_argument("invalid orbital index: " + s);
    std::size_t n = 0;
    const auto value = std::stoull(s, &n);
    if (n != s.size() || value == 0 || value > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("invalid orbital index: " + s);
    return static_cast<std::size_t>(value);
}

// Resolve an absolute MO index or a HOMO/LUMO-relative index (1-based).
std::size_t endpoint(const std::string& s, const Molecule& mol, Spin channel) {
    if (s.empty()) throw std::invalid_argument("empty orbital selector");
    if (std::isdigit(static_cast<unsigned char>(s[0]))) return number(s);
    std::size_t base = 0, offset = 0;
    if (s.rfind("homo", 0) == 0) {
        base = static_cast<std::size_t>(mol.occupied_count(channel));
        offset = 4;
    } else if (s.rfind("lumo", 0) == 0) {
        base = static_cast<std::size_t>(mol.occupied_count(channel)) + 1;
        offset = 4;
    } else
        throw std::invalid_argument("unknown MO selector: " + s);
    if (offset == s.size()) return base;
    const char sign = s[offset];
    if (sign != '+' && sign != '-') throw std::invalid_argument("invalid MO offset: " + s);
    const auto step = number(s.substr(offset + 1));
    if (sign == '-' && step >= base) throw std::out_of_range("MO offset below first orbital: " + s);
    if (sign == '+') {
        if (step > std::numeric_limits<std::size_t>::max() - base)
            throw std::out_of_range("MO index overflow");
        return base + step;
    }
    return base - step;
}

// Expand inclusive ranges and de-duplicate MO selections in user order.
std::vector<std::size_t> select_mos(const std::string& selector, const Molecule& mol,
                                    Spin channel) {
    std::vector<std::size_t> result;
    std::set<std::size_t> seen;
    std::string expanded = selector == "frontier" ? "homo,lumo" : selector;
    std::size_t start = 0;
    while (start < expanded.size()) {
        const auto comma = expanded.find(',', start);
        const auto token =
            expanded.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        const auto colon = token.find(':');
        auto first = endpoint(token.substr(0, colon), mol, channel);
        auto last =
            colon == std::string::npos ? first : endpoint(token.substr(colon + 1), mol, channel);
        if (first > last || last > mol.orbital_count(channel))
            throw std::out_of_range("MO selector outside available orbitals: " + token);
        if (last - first > 10000) throw std::invalid_argument("MO range too large");
        for (auto n = first; n <= last; ++n)
            if (seen.insert(n).second) result.push_back(n);
        if (comma == std::string::npos) break;
        start = comma + 1;
        if (start == expanded.size()) throw std::invalid_argument("trailing comma in MO selector");
    }
    if (result.empty()) throw std::invalid_argument("empty MO selector");
    return result;
}

std::string spin_name(Spin s) {
    return s == Spin::Alpha ? "alpha" : "beta";
}

void list_mos(const Molecule& mol, Spin channel) {
    const auto count = mol.orbital_count(channel);
    const auto& energies = mol.orbital_energies(channel);
    const auto occupied = mol.occupied_count(channel);
    const bool restricted = !mol.unrestricted();
    std::cout << "Spin: " << (restricted ? "restricted" : spin_name(channel))
              << "  orbitals: " << count << "  occupied: " << occupied << "\n";
    std::cout << " MO   Occupation   Energy (hartree)   Label\n";
    for (std::size_t i = 1; i <= count; ++i) {
        std::string label;
        if (i == static_cast<std::size_t>(occupied)) label = "HOMO";
        if (i == static_cast<std::size_t>(occupied) + 1) label = "LUMO";
        std::cout << std::setw(4) << i << "   "
                  << (i <= static_cast<std::size_t>(occupied) ? (restricted ? 2 : 1) : 0)
                  << "            ";
        if (i <= energies.size())
            std::cout << std::setw(15) << std::fixed << std::setprecision(8) << energies[i - 1];
        else
            std::cout << std::setw(15) << "n/a";
        std::cout << "   " << label << "\n";
    }
}

void write_orbitals(const Molecule& mol, const chem::kumo::Grid& grid, const Options& o,
                    Spin channel) {
    const auto indices = select_mos(o.selector, mol, channel);
    auto fields = mol.orbitals(grid, indices, channel);
    for (std::size_t j = 0; j < indices.size(); ++j) {
        const auto n = indices[j];
        const auto filename =
            o.prefix + "_" + spin_name(channel) + "_mo_" + std::to_string(n) + ".cube";
        fields[j].write_orbital_cube(filename, n,
                                     "Molecular orbital amplitude (signed), " + spin_name(channel) +
                                         " MO " + std::to_string(n));
        std::cout << filename << "\n";
    }
}

// Density property requested with --property, with its Cube filename suffix.
struct DensityProperty {
    PropertyKind kind;
    std::string name;
    std::string description; // Cube comment line, with units
};

// Parse --property: a comma-separated list, duplicates ignored, user order kept.
std::vector<DensityProperty> select_properties(const std::string& list) {
    const DensityProperty known[] = {
        {PropertyKind::Density, "density", "Electron density (e/bohr^3)"},
        {PropertyKind::Gradient, "gradient",
         "Electron density gradient norm |grad rho| (e/bohr^4)"},
        {PropertyKind::Laplacian, "laplacian", "Electron density Laplacian (e/bohr^5)"},
        {PropertyKind::Tau, "tau", "Positive kinetic energy density tau (hartree/bohr^3)"}};
    std::vector<DensityProperty> result;
    std::size_t start = 0;
    while (true) {
        const auto comma = list.find(',', start);
        const auto token =
            list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        const auto* match = std::find_if(std::begin(known), std::end(known),
                                         [&](const DensityProperty& p) { return p.name == token; });
        if (match == std::end(known))
            throw std::invalid_argument("unknown --property '" + token +
                                        "' (use density, gradient, laplacian or tau)");
        if (std::none_of(result.begin(), result.end(),
                         [&](const DensityProperty& p) { return p.kind == match->kind; }))
            result.push_back(*match);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return result;
}

std::string component_label(Spin component) {
    switch (component) {
    case Spin::Total: return "total";
    case Spin::Alpha: return "alpha";
    case Spin::Beta: return "beta";
    default: return "alpha-beta difference";
    }
}

chem::kumo::DensityDescriptor descriptor(Spin component) {
    switch (component) {
    case Spin::Total: return chem::kumo::density::total;
    case Spin::Alpha: return chem::kumo::density::alpha;
    case Spin::Beta: return chem::kumo::density::beta;
    default: return chem::kumo::density::spin;
    }
}

chem::kumo::Property density_property(Spin component, PropertyKind kind) {
    const auto d = descriptor(component);
    switch (kind) {
    case PropertyKind::Gradient: return d.gradient();
    case PropertyKind::Laplacian: return d.laplacian();
    case PropertyKind::Tau: return d.tau();
    default: return d.value();
    }
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

std::string backend_label(chem::kumo::Backend backend) {
    return backend == chem::kumo::Backend::AO_MO ? "AO-MO contraction"
                                                 : "AO density-matrix contraction";
}

// Evaluate every requested property of every density component in one
// Evaluator request and write one Cube file each: first the densities
// (PREFIX_<component>.cube), then the derivatives
// (PREFIX_<component>_<property>.cube). The gradient is written as its norm.
// Each Cube comment names the quantity, its units, the component and the
// contraction actually used.
void write_density_properties(const Molecule& mol, const chem::kumo::Grid& grid,
                              const std::vector<std::pair<Spin, std::string>>& components,
                              const std::vector<DensityProperty>& properties,
                              chem::kumo::Backend backend, const std::string& prefix) {
    chem::kumo::Request request;
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            request.add(density_property(component, property.kind));
    chem::kumo::EvaluationOptions options;
    options.backend = backend;
    const auto results = chem::kumo::Evaluator(mol, options).evaluate(grid, request);
    const auto used = backend_label(results.backend_used());
    const auto write = [&](Spin component, const std::string& suffix,
                           const DensityProperty& property) {
        const auto& field = results.get(density_property(component, property.kind));
        std::vector<double> values(field.points);
        if (property.kind == PropertyKind::Gradient) {
            for (std::size_t q = 0; q < field.points; ++q)
                values[q] = std::sqrt(field.component(0)[q] * field.component(0)[q] +
                                      field.component(1)[q] * field.component(1)[q] +
                                      field.component(2)[q] * field.component(2)[q]);
        } else {
            const auto source = field.component();
            std::copy(source.begin(), source.end(), values.begin());
        }
        const auto filename = prefix + "_" + suffix +
                              (property.kind == PropertyKind::Density ? "" : "_" + property.name) +
                              ".cube";
        chem::kumo::Field(grid, mol.atoms(), std::move(values))
            .write_cube(filename,
                        property.description + ", " + component_label(component) + ", " + used);
        std::cout << filename << "\n";
    };
    for (const auto& property : properties)
        if (property.kind == PropertyKind::Density)
            for (const auto& [component, suffix] : components)
                write(component, suffix, property);
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            if (property.kind != PropertyKind::Density) write(component, suffix, property);
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") {
            usage();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--version") {
            std::cout << "kumo " CHEM_KUMO_VERSION "\n";
            return 0;
        }
        const auto o = parse(argc, argv);
        const auto mol = Molecule::from_fchk(o.input);
        if (o.list) {
            if (o.spin.empty() || o.spin == "alpha" || o.spin == "both") list_mos(mol, Spin::Alpha);
            if (o.spin == "beta" || o.spin == "both") list_mos(mol, Spin::Beta);
            if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
                throw std::invalid_argument("--list-mo accepts alpha, beta, or both");
            return 0;
        }
        const auto grid = chem::kumo::Grid::around(mol.atoms(), o.margin, o.spacing);
        if (o.mo_mode) {
            if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
                throw std::invalid_argument("--mo accepts --spin alpha, beta, or both");
            if (o.spin.empty() || o.spin == "alpha" || o.spin == "both")
                write_orbitals(mol, grid, o, Spin::Alpha);
            if (o.spin == "beta" || o.spin == "both") write_orbitals(mol, grid, o, Spin::Beta);
        } else {
            const auto spin = o.spin.empty() ? "total" : o.spin;
            std::vector<std::pair<Spin, std::string>> components;
            if (spin == "total" || spin == "all") components.emplace_back(Spin::Total, "density");
            if (spin == "alpha" || spin == "all") components.emplace_back(Spin::Alpha, "alpha");
            if (spin == "beta" || spin == "all") components.emplace_back(Spin::Beta, "beta");
            if (spin == "difference" || spin == "all")
                components.emplace_back(Spin::Difference, "spin");
            if (components.empty())
                throw std::invalid_argument(
                    "density --spin must be total, alpha, beta, difference, or all");
            const auto properties =
                select_properties(o.properties.empty() ? "density" : o.properties);
            write_density_properties(mol, grid, components, properties,
                                     select_backend(o.backend, mol), o.prefix);
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "kumo: " << e.what() << "\nTry 'kumo --help'.\n";
        return 1;
    }
}
