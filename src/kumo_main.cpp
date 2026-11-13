#include "chem/kumo.hpp"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <io.h>
#define KUMO_ISATTY(stream) _isatty(_fileno(stream))
#else
#include <unistd.h>
#define KUMO_ISATTY(stream) isatty(fileno(stream))
#endif

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
    double esp_tolerance = 1e-8; // Evaluator default
    bool esp_tolerance_set = false;
    double isovalue = 1e-3; // --esp-surface and --alie-surface density isosurface
    bool isovalue_set = false;
    bool spacing_set = false;
    bool margin_set = false;
    bool esp_surface = false;
    bool alie_surface = false;
    bool mo_mode = false;
    bool list = false;
    bool verbose = false;
    bool output_set = false;
};

void usage() {
    std::cout << "kumo " CHEM_KUMO_VERSION
                 " - ground-state electron densities and molecular orbitals\n\n"
              << "Usage: kumo FILE.fchk[.gz] [options]\n\n"
              << "  (default)                  Write total electron density\n"
              << "  --spin total|alpha|beta|difference|all   Density component(s)\n"
              << "  --property LIST            Density properties to write, comma-separated\n"
              << "                             (default density):\n"
              << "                             density, gradient (the norm |grad rho|),\n"
              << "                             laplacian, tau, kinetic-k (any --spin);\n"
              << "                             elf, lol (total, alpha or beta);\n"
              << "                             rdg, sign-lambda2-rho, iri, dori,\n"
              << "                             potential-energy, energy-density, shannon,\n"
              << "                             fisher, esp (the electrostatic potential)\n"
              << "                             (total only); spin-polarization (difference)\n"
              << "  --esp-tolerance HARTREE    Error bound of the electrostatic potential at\n"
              << "                             each point (default 1e-8; 0 is exact and slowest)\n"
              << "  --esp-surface              Print the electrostatic potential statistics on\n"
              << "                             a density isosurface (no Cube files; grid\n"
              << "                             spacing default 0.25 bohr)\n"
              << "  --alie-surface             Print the average local ionization energy\n"
              << "                             statistics on a density isosurface (as\n"
              << "                             --esp-surface)\n"
              << "  --isovalue RHO             Isosurface of --esp-surface or --alie-surface\n"
              << "                             (default 0.001)\n"
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
              << "  -v, --verbose              Report progress on standard error\n"
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

// Progress messages for --verbose, written to standard error so that standard
// output keeps listing only the files written. On a terminal, evaluations
// also show a progress bar that is redrawn in place; otherwise only the
// start and the end of each step are reported.
class Console {
public:
    using clock = std::chrono::steady_clock;

    explicit Console(bool verbose)
        : verbose_(verbose), terminal_(verbose && KUMO_ISATTY(stderr)), start_(clock::now()) {}

    bool verbose() const noexcept { return verbose_; }

    void say(const std::string& line) const {
        if (verbose_) std::cerr << line << '\n' << std::flush;
    }

    /// Announce a step and start timing it.
    void begin(const std::string& step) {
        say(step);
        step_start_ = clock::now();
        last_percent_ = -1;
        bar_drawn_ = false;
    }

    /// Progress callback for the current step (empty unless on a terminal).
    chem::kumo::ProgressCallback progress() {
        if (!terminal_) return {};
        return [this](std::size_t done, std::size_t total) { draw(done, total); };
    }

    /// End the current step with a summary such as "Evaluated in".
    void end(const std::string& summary) {
        if (bar_drawn_) std::cerr << '\n';
        say("  " + summary + " " + seconds(since(step_start_)));
    }

    void finish() const { say("Finished in " + seconds(since(start_))); }

    static std::string seconds(double s) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(s < 10 ? 2 : 1) << s << " s";
        return out.str();
    }

    /// 1543260 -> "1,543,260".
    static std::string grouped(std::size_t n) {
        auto digits = std::to_string(n);
        for (auto i = static_cast<std::ptrdiff_t>(digits.size()) - 3; i > 0; i -= 3)
            digits.insert(static_cast<std::size_t>(i), ",");
        return digits;
    }

private:
    static double since(clock::time_point t) {
        return std::chrono::duration<double>(clock::now() - t).count();
    }

    // "  [==========>         ]  52%  5.2 s elapsed, about 4.8 s left"
    void draw(std::size_t done, std::size_t total) {
        if (total == 0) return;
        const int percent = static_cast<int>(100 * done / total);
        if (percent == last_percent_) return;
        last_percent_ = percent;
        constexpr int width = 30;
        const int filled = static_cast<int>(width * done / total);
        std::string bar(static_cast<std::size_t>(filled), '=');
        if (filled < width) bar += '>';
        bar.resize(width, ' ');
        const double elapsed = since(step_start_);
        std::ostringstream line;
        line << "\r  [" << bar << "] " << std::setw(3) << percent << "%  " << seconds(elapsed)
             << " elapsed";
        if (done > 0 && done < total)
            line << ", about "
                 << seconds(elapsed * static_cast<double>(total - done) / static_cast<double>(done))
                 << " left";
        line << "    ";
        std::cerr << line.str() << std::flush;
        bar_drawn_ = true;
    }

    bool verbose_;
    bool terminal_;
    clock::time_point start_;
    clock::time_point step_start_ = clock::now();
    int last_percent_ = -1;
    bool bar_drawn_ = false;
};

// "a", "a and b", "a, b and c".
std::string join(const std::vector<std::string>& items) {
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i > 0) out += i + 1 == items.size() ? " and " : ", ";
        out += items[i];
    }
    return out;
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
        } else if (a == "--esp-tolerance") {
            if (o.esp_tolerance_set) throw std::invalid_argument("--esp-tolerance specified twice");
            o.esp_tolerance = real(value());
            o.esp_tolerance_set = true;
        } else if (a == "--esp-surface")
            o.esp_surface = true;
        else if (a == "--alie-surface")
            o.alie_surface = true;
        else if (a == "--isovalue") {
            if (o.isovalue_set) throw std::invalid_argument("--isovalue specified twice");
            o.isovalue = real(value());
            o.isovalue_set = true;
        } else if (a == "--spacing") {
            o.spacing = real(value());
            o.spacing_set = true;
        } else if (a == "--margin") {
            o.margin = real(value());
            o.margin_set = true;
        } else if (a == "--output") {
            o.prefix = value();
            o.output_set = true;
        } else if (a == "--list-mo")
            o.list = true;
        else if (a == "-v" || a == "--verbose")
            o.verbose = true;
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
    if (o.esp_surface && o.alie_surface)
        throw std::invalid_argument("choose --esp-surface or --alie-surface");
    if (o.esp_surface || o.alie_surface) {
        const std::string option = o.esp_surface ? "--esp-surface" : "--alie-surface";
        if (o.mo_mode || o.list || !o.properties.empty() || !o.spin.empty() || !o.backend.empty() ||
            o.margin_set || o.output_set)
            throw std::invalid_argument(option + " does not combine with --mo, --list-mo, "
                                                 "--property, --spin, --backend, --margin or "
                                                 "--output");
        if (o.alie_surface && o.esp_tolerance_set)
            throw std::invalid_argument("--esp-tolerance does not apply to --alie-surface");
        if (!o.spacing_set) o.spacing = 0.25;
    } else if (o.isovalue_set)
        throw std::invalid_argument("--isovalue applies only to --esp-surface or --alie-surface");
    if (!(o.spacing > 0) || o.margin < 0)
        throw std::invalid_argument("invalid grid spacing or margin");
    if (!std::isfinite(o.esp_tolerance) || o.esp_tolerance < 0)
        throw std::invalid_argument("--esp-tolerance must be finite and nonnegative");
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
                    Spin channel, Console& console) {
    const auto indices = select_mos(o.selector, mol, channel);
    std::vector<std::string> numbers;
    for (auto n : indices)
        numbers.push_back(std::to_string(n));
    console.begin("Evaluating " + spin_name(channel) + " molecular orbital" +
                  (indices.size() > 1 ? "s " : " ") + join(numbers));
    auto fields = mol.orbitals(grid, indices, channel, console.progress());
    console.end("Evaluated in");
    console.begin("Writing " + std::to_string(indices.size()) + " Cube file" +
                  (indices.size() > 1 ? "s" : ""));
    for (std::size_t j = 0; j < indices.size(); ++j) {
        const auto n = indices[j];
        const auto filename =
            o.prefix + "_" + spin_name(channel) + "_mo_" + std::to_string(n) + ".cube";
        fields[j].write_orbital_cube(filename, n,
                                     "Molecular orbital amplitude (signed), " + spin_name(channel) +
                                         " MO " + std::to_string(n));
        std::cout << filename << "\n";
    }
    console.end("Written in");
}

// Density property requested with --property, with its Cube filename suffix.
struct DensityProperty {
    PropertyKind kind;
    std::string name;
    std::string description; // Cube comment line, with units
    std::string summary;     // in "Evaluating the <summary> of the total electron density"
};

// Parse --property: a comma-separated list, duplicates ignored, user order kept.
std::vector<DensityProperty> select_properties(const std::string& list) {
    const DensityProperty known[] = {
        {PropertyKind::Density, "density", "Electron density (e/bohr^3)", "value"},
        {PropertyKind::Gradient, "gradient", "Electron density gradient norm |grad rho| (e/bohr^4)",
         "gradient norm"},
        {PropertyKind::Laplacian, "laplacian", "Electron density Laplacian (e/bohr^5)",
         "Laplacian"},
        {PropertyKind::Tau, "tau", "Positive kinetic energy density tau (hartree/bohr^3)",
         "kinetic energy density tau"},
        {PropertyKind::Elf, "elf", "Electron localization function ELF (dimensionless)", "ELF"},
        {PropertyKind::Lol, "lol", "Localized orbital locator LOL (dimensionless)", "LOL"},
        {PropertyKind::ReducedGradient, "rdg",
         "Reduced density gradient s (dimensionless, capped at 100)", "reduced density gradient"},
        {PropertyKind::SignLambda2Density, "sign-lambda2-rho",
         "sign(lambda2) * electron density (e/bohr^3)", "sign(lambda2) rho"},
        {PropertyKind::InteractionRegionIndicator, "iri",
         "Interaction region indicator |grad rho| / rho^1.1 (capped at 100)",
         "interaction region indicator"},
        {PropertyKind::DensityOverlapRegionsIndicator, "dori",
         "Density overlap regions indicator DORI (dimensionless)",
         "density overlap regions indicator"},
        {PropertyKind::HamiltonianKineticEnergy, "kinetic-k",
         "Hamiltonian kinetic energy density K = tau - lap rho / 4 (hartree/bohr^3)",
         "kinetic energy density K"},
        {PropertyKind::PotentialEnergyDensity, "potential-energy",
         "Potential energy density V = lap rho / 4 - 2 tau (hartree/bohr^3)",
         "potential energy density"},
        {PropertyKind::EnergyDensity, "energy-density",
         "Electronic energy density H = lap rho / 4 - tau (hartree/bohr^3)", "energy density"},
        {PropertyKind::SpinPolarization, "spin-polarization",
         "Spin polarization (rho_alpha - rho_beta) / rho (dimensionless)", "spin polarization"},
        {PropertyKind::ShannonEntropyDensity, "shannon",
         "Shannon entropy density -rho ln rho (e/bohr^3)", "Shannon entropy density"},
        {PropertyKind::FisherInformationDensity, "fisher",
         "Fisher information density |grad rho|^2 / rho (e/bohr^5)", "Fisher information density"},
        {PropertyKind::ElectrostaticPotential, "esp",
         "Molecular electrostatic potential, nuclear + electronic (hartree/e)",
         "electrostatic potential"}};
    std::vector<DensityProperty> result;
    std::size_t start = 0;
    while (true) {
        const auto comma = list.find(',', start);
        const auto token =
            list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        const auto* match = std::find_if(std::begin(known), std::end(known),
                                         [&](const DensityProperty& p) { return p.name == token; });
        if (match == std::end(known)) {
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
    case PropertyKind::Elf: return d.elf();
    case PropertyKind::ReducedGradient: return d.reduced_gradient();
    case PropertyKind::SignLambda2Density: return d.sign_lambda2_density();
    case PropertyKind::ElectrostaticPotential: return d.electrostatic_potential();
    case PropertyKind::Lol: return d.lol();
    case PropertyKind::InteractionRegionIndicator: return d.interaction_region_indicator();
    case PropertyKind::DensityOverlapRegionsIndicator: return d.density_overlap_regions_indicator();
    case PropertyKind::HamiltonianKineticEnergy: return d.hamiltonian_kinetic_energy();
    case PropertyKind::PotentialEnergyDensity: return d.potential_energy_density();
    case PropertyKind::EnergyDensity: return d.energy_density();
    case PropertyKind::SpinPolarization: return d.spin_polarization();
    case PropertyKind::ShannonEntropyDensity: return d.shannon_entropy_density();
    case PropertyKind::FisherInformationDensity: return d.fisher_information_density();
    default: return d.value();
    }
}

// ELF is not defined for the spin difference, and the reduced gradient,
// sign(lambda2) rho and the electrostatic potential only for the total
// density.
void check_components(const std::vector<std::pair<Spin, std::string>>& components,
                      const std::vector<DensityProperty>& properties) {
    for (const auto& property : properties)
        for (const auto& component : components) {
            const Spin spin = component.first;
            if ((property.kind == PropertyKind::Elf || property.kind == PropertyKind::Lol) &&
                spin == Spin::Difference)
                throw std::invalid_argument("--property " + property.name +
                                            " is not defined for --spin difference (or all); "
                                            "use total, alpha or beta");
            if (property.kind == PropertyKind::SpinPolarization && spin != Spin::Difference)
                throw std::invalid_argument(
                    "--property spin-polarization is defined only for --spin difference");
            const bool total_only = property.kind == PropertyKind::ReducedGradient ||
                                    property.kind == PropertyKind::SignLambda2Density ||
                                    property.kind == PropertyKind::ElectrostaticPotential ||
                                    property.kind == PropertyKind::InteractionRegionIndicator ||
                                    property.kind == PropertyKind::DensityOverlapRegionsIndicator ||
                                    property.kind == PropertyKind::PotentialEnergyDensity ||
                                    property.kind == PropertyKind::EnergyDensity ||
                                    property.kind == PropertyKind::ShannonEntropyDensity ||
                                    property.kind == PropertyKind::FisherInformationDensity;
            if (total_only && spin != Spin::Total)
                throw std::invalid_argument("--property " + property.name +
                                            " is defined only for --spin total");
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

// "Evaluating the total electron density", or "Evaluating the gradient norm
// and Laplacian of the total and alpha electron densities".
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

// Evaluate every requested property of every density component in one
// Evaluator request and write one Cube file each: first the densities
// (PREFIX_<component>.cube), then the derivatives
// (PREFIX_<component>_<property>.cube). The gradient is written as its norm.
// Each Cube comment names the quantity, its units, the component and the
// contraction actually used.
void write_density_properties(const Molecule& mol, const chem::kumo::Grid& grid,
                              const std::vector<std::pair<Spin, std::string>>& components,
                              const std::vector<DensityProperty>& properties,
                              chem::kumo::Backend backend, double esp_tolerance,
                              const std::string& prefix, Console& console) {
    chem::kumo::Request request;
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            request.add(density_property(component, property.kind));
    console.begin(evaluation_summary(components, properties));
    chem::kumo::EvaluationOptions options;
    options.backend = backend;
    options.potential_tolerance = esp_tolerance;
    options.progress = console.progress();
    const auto results = chem::kumo::Evaluator(mol, options).evaluate(grid, request);
    console.end("Evaluated in");
    const auto files = components.size() * properties.size();
    console.begin("Writing " + std::to_string(files) + " Cube file" + (files > 1 ? "s" : ""));
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
        // The electrostatic potential always uses the AO density matrix.
        std::ostringstream bound;
        bound << "error <= " << esp_tolerance << " hartree/e";
        const auto method = property.kind != PropertyKind::ElectrostaticPotential ? used
                            : esp_tolerance == 0.0
                                ? std::string("exact Coulomb integrals over the AO density matrix")
                                : "Coulomb integrals over the AO density matrix, " + bound.str();
        chem::kumo::Field(grid, mol.atoms(), std::move(values))
            .write_cube(filename,
                        property.description + ", " + component_label(component) + ", " + method);
        std::cout << filename << "\n";
    };
    for (const auto& property : properties)
        if (property.kind == PropertyKind::Density)
            for (const auto& [component, suffix] : components)
                write(component, suffix, property);
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            if (property.kind != PropertyKind::Density) write(component, suffix, property);
    console.end("Written in");
}

// "60 atoms, 1410 basis functions; restricted closed-shell SCF, 90 doubly
// occupied orbitals".
std::string molecule_summary(const Molecule& mol) {
    const auto na = mol.occupied_count(Spin::Alpha), nb = mol.occupied_count(Spin::Beta);
    std::string reference;
    if (mol.unrestricted())
        reference = "unrestricted SCF, " + std::to_string(na) + " alpha and " + std::to_string(nb) +
                    " beta occupied orbitals";
    else if (na == nb)
        reference =
            "restricted closed-shell SCF, " + std::to_string(na) + " doubly occupied orbitals";
    else
        reference = "restricted open-shell SCF, " + std::to_string(na) + " alpha and " +
                    std::to_string(nb) + " beta occupied orbitals";
    return std::to_string(mol.atoms().size()) + " atoms, " + std::to_string(mol.basis_size()) +
           " basis functions; " + reference;
}

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

// --esp-surface: statistics of the electrostatic potential on the density
// isosurface, as a report on standard output.
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
}

// --alie-surface: statistics of the average local ionization energy on the
// density isosurface, as a report on standard output.
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
}

// "Grid: 170 x 178 x 51 = 1,543,260 points, spacing 0.200 bohr, margin 5.000 bohr".
std::string grid_summary(const chem::kumo::Grid& grid, const Options& o) {
    std::ostringstream out;
    out << "Grid: " << grid.counts[0] << " x " << grid.counts[1] << " x " << grid.counts[2] << " = "
        << Console::grouped(static_cast<std::size_t>(grid.counts[0]) *
                            static_cast<std::size_t>(grid.counts[1]) *
                            static_cast<std::size_t>(grid.counts[2]))
        << " points, spacing " << std::fixed << std::setprecision(3) << o.spacing
        << " bohr, margin " << o.margin << " bohr";
    return out.str();
}

// Result of the load-time check of the occupied orbitals against the FCHK
// density matrices, which decides whether AO-MO can be used.
std::string orbital_check_summary(const Molecule& mol) {
    const auto& info = mol.electronic_structure_info();
    if (info.orbital_status == chem::kumo::OrbitalRepresentationStatus::Unavailable)
        return "Occupied orbitals: not available in the FCHK, so AO-MO cannot be used";
    std::ostringstream out;
    out << std::setprecision(2)
        << "Occupied orbitals vs FCHK density matrices: error relative to the total matrix "
        << info.total_relative_error << " (total), " << info.spin_relative_error
        << " (spin); AO-MO needs both <= " << chem::kumo::orbital_check_relative_tolerance;
    if (!info.ao_mo_eligible()) out << " (not met)";
    return out.str();
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
        Console console(o.verbose);
        console.say("kumo " CHEM_KUMO_VERSION);
        console.begin("Reading formatted checkpoint " + o.input.string());
        const auto mol = Molecule::from_fchk(o.input);
        console.end("Read in");
        console.say("  " + molecule_summary(mol));
        console.say("  " + orbital_check_summary(mol));
        if (o.list) {
            if (o.spin.empty() || o.spin == "alpha" || o.spin == "both") list_mos(mol, Spin::Alpha);
            if (o.spin == "beta" || o.spin == "both") list_mos(mol, Spin::Beta);
            if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
                throw std::invalid_argument("--list-mo accepts alpha, beta, or both");
            return 0;
        }
        if (o.esp_surface || o.alie_surface) {
            if (o.esp_surface)
                report_esp_surface(mol, o, console);
            else
                report_alie_surface(mol, o, console);
            console.finish();
            return 0;
        }
        const auto grid = chem::kumo::Grid::around(mol.atoms(), o.margin, o.spacing);
        console.say(grid_summary(grid, o));
        if (o.mo_mode) {
            if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
                throw std::invalid_argument("--mo accepts --spin alpha, beta, or both");
            if (o.spin.empty() || o.spin == "alpha" || o.spin == "both")
                write_orbitals(mol, grid, o, Spin::Alpha, console);
            if (o.spin == "beta" || o.spin == "both")
                write_orbitals(mol, grid, o, Spin::Beta, console);
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
            check_components(components, properties);
            if (o.esp_tolerance_set &&
                std::none_of(properties.begin(), properties.end(), [](const DensityProperty& p) {
                    return p.kind == PropertyKind::ElectrostaticPotential;
                }))
                throw std::invalid_argument("--esp-tolerance applies only to --property esp");
            const auto backend = select_backend(o.backend, mol);
            console.say(backend_summary(backend, mol));
            write_density_properties(mol, grid, components, properties, backend, o.esp_tolerance,
                                     o.prefix, console);
        }
        console.finish();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "kumo: " << e.what() << "\nTry 'kumo --help'.\n";
        return 1;
    }
}
