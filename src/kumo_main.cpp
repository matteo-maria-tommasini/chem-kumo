#include "chem/kumo.hpp"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
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
    bool surface_map = false; // --surface-map
    bool ply = false;         // --ply
    bool mo_mode = false;
    bool list = false;
    bool verbose = false;
    bool output_set = false;
    std::vector<std::string> points; // --point, repeatable
    std::string line;                // --line A:B
    std::string plane;               // --plane A:B:C
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
              << "                             fisher, esp (the electrostatic potential),\n"
              << "                             alie, lea, leae, local-electronegativity,\n"
              << "                             local-hardness (orbital-energy fields of SCF\n"
              << "                             orbitals) (total only); spin-polarization\n"
              << "                             (difference)\n"
              << "  --esp-tolerance HARTREE    Error bound of the electrostatic potential at\n"
              << "                             each point (default 1e-8; 0 is exact and slowest)\n"
              << "  --esp-surface              Print the electrostatic potential statistics on\n"
              << "                             a density isosurface (no Cube files; grid\n"
              << "                             spacing default 0.25 bohr)\n"
              << "  --alie-surface             Print the average local ionization energy\n"
              << "                             statistics on a density isosurface (as\n"
              << "                             --esp-surface)\n"
              << "  --surface-map              Print the statistics of each --property (with\n"
              << "                             --spin) on a density isosurface (as\n"
              << "                             --esp-surface)\n"
              << "  --isovalue RHO             Isosurface of --esp-surface, --alie-surface or\n"
              << "                             --surface-map (default 0.001)\n"
              << "  --ply                      With a surface option, also write the surface\n"
              << "                             and its values to PREFIX_surface.ply (bohr)\n"
              << "  --point X,Y,Z|@N           Print the density properties at a point, in bohr\n"
              << "                             or at atom N (1-based); repeatable\n"
              << "  --line A:B                 Write the properties along the line through A\n"
              << "                             and B (X,Y,Z or @N) to PREFIX_line.dat; spacing\n"
              << "                             default 0.01 bohr, margin (extension beyond A\n"
              << "                             and B) default 0\n"
              << "  --plane A:B:C              Write the properties on the plane of A, B and C\n"
              << "                             to PREFIX_plane.dat; spacing default 0.10 bohr,\n"
              << "                             margin around the three points default 5.0 bohr\n"
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
        else if (a == "--surface-map")
            o.surface_map = true;
        else if (a == "--ply")
            o.ply = true;
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
        } else if (a == "--point")
            o.points.push_back(value());
        else if (a == "--line") {
            if (!o.line.empty()) throw std::invalid_argument("--line specified twice");
            o.line = value();
            if (o.line.empty()) throw std::invalid_argument("empty --line");
        } else if (a == "--plane") {
            if (!o.plane.empty()) throw std::invalid_argument("--plane specified twice");
            o.plane = value();
            if (o.plane.empty()) throw std::invalid_argument("empty --plane");
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
    const int samplings =
        (o.points.empty() ? 0 : 1) + (o.line.empty() ? 0 : 1) + (o.plane.empty() ? 0 : 1);
    if (samplings > 1) throw std::invalid_argument("choose one of --point, --line or --plane");
    if (samplings == 1) {
        if (o.mo_mode || o.list || o.esp_surface || o.alie_surface)
            throw std::invalid_argument("--point, --line and --plane apply only to density "
                                        "properties, not to --mo, --list-mo, --esp-surface or "
                                        "--alie-surface");
        if (!o.points.empty() && (o.spacing_set || o.margin_set || o.output_set))
            throw std::invalid_argument(
                "--point prints to standard output and does not use --spacing, --margin or "
                "--output");
        if (!o.line.empty()) {
            if (!o.spacing_set) o.spacing = 0.01;
            if (!o.margin_set) o.margin = 0.0;
        }
        if (!o.plane.empty() && !o.spacing_set) o.spacing = 0.10;
    }
    if ((o.esp_surface ? 1 : 0) + (o.alie_surface ? 1 : 0) + (o.surface_map ? 1 : 0) > 1)
        throw std::invalid_argument("choose one of --esp-surface, --alie-surface or --surface-map");
    if (o.surface_map) {
        if (samplings > 0 || o.mo_mode || o.list || o.margin_set)
            throw std::invalid_argument("--surface-map does not combine with --point, --line, "
                                        "--plane, --mo, --list-mo or --margin");
        if (o.properties.empty())
            throw std::invalid_argument("--surface-map needs --property LIST");
        if (o.output_set && !o.ply)
            throw std::invalid_argument("--output applies to --surface-map only with --ply");
        if (!o.spacing_set) o.spacing = 0.25;
    } else if (o.esp_surface || o.alie_surface) {
        const std::string option = o.esp_surface ? "--esp-surface" : "--alie-surface";
        if (o.mo_mode || o.list || !o.properties.empty() || !o.spin.empty() || !o.backend.empty() ||
            o.margin_set || (o.output_set && !o.ply))
            throw std::invalid_argument(option + " does not combine with --mo, --list-mo, "
                                                 "--property, --spin, --backend, --margin, or "
                                                 "--output without --ply");
        if (o.alie_surface && o.esp_tolerance_set)
            throw std::invalid_argument("--esp-tolerance does not apply to --alie-surface");
        if (!o.spacing_set) o.spacing = 0.25;
    } else if (o.isovalue_set)
        throw std::invalid_argument(
            "--isovalue applies only to --esp-surface, --alie-surface or --surface-map");
    else if (o.ply)
        throw std::invalid_argument(
            "--ply applies only to --esp-surface, --alie-surface or --surface-map");
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
         "electrostatic potential"},
        {PropertyKind::AverageLocalIonizationEnergy, "alie",
         "Average local ionization energy (hartree)", "average local ionization energy"},
        {PropertyKind::LocalElectronAffinity, "lea", "Local electron affinity (hartree)",
         "local electron affinity"},
        {PropertyKind::LocalElectronAttachmentEnergy, "leae",
         "Local electron attachment energy (hartree)", "local electron attachment energy"},
        {PropertyKind::LocalElectronegativity, "local-electronegativity",
         "Local electronegativity (ALIE + LEA) / 2 (hartree)", "local electronegativity"},
        {PropertyKind::LocalHardness, "local-hardness", "Local hardness (ALIE - LEA) / 2 (hartree)",
         "local hardness"}};
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
    case PropertyKind::AverageLocalIonizationEnergy: return d.average_local_ionization_energy();
    case PropertyKind::LocalElectronAffinity: return d.local_electron_affinity();
    case PropertyKind::LocalElectronAttachmentEnergy: return d.local_electron_attachment_energy();
    case PropertyKind::LocalElectronegativity: return d.local_electronegativity();
    case PropertyKind::LocalHardness: return d.local_hardness();
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
                                    property.kind == PropertyKind::FisherInformationDensity ||
                                    property.kind == PropertyKind::AverageLocalIonizationEnergy ||
                                    property.kind == PropertyKind::LocalElectronAffinity ||
                                    property.kind == PropertyKind::LocalElectronAttachmentEnergy ||
                                    property.kind == PropertyKind::LocalElectronegativity ||
                                    property.kind == PropertyKind::LocalHardness;
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

// Values of a property at every point, as written to Cube files and tables:
// the gradient as its norm, the other properties as they are.
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

// How a property was computed, for Cube comments and table headers. The
// electrostatic potential always uses the AO density matrix.
std::string method_label(PropertyKind kind, chem::kumo::Backend used, double esp_tolerance) {
    if (kind == PropertyKind::AverageLocalIonizationEnergy ||
        kind == PropertyKind::LocalElectronAffinity ||
        kind == PropertyKind::LocalElectronAttachmentEnergy ||
        kind == PropertyKind::LocalElectronegativity || kind == PropertyKind::LocalHardness)
        return "SCF orbitals and orbital energies";
    if (kind != PropertyKind::ElectrostaticPotential) return backend_label(used);
    if (esp_tolerance == 0.0) return "exact Coulomb integrals over the AO density matrix";
    std::ostringstream bound;
    bound << "error <= " << esp_tolerance << " hartree/e";
    return "Coulomb integrals over the AO density matrix, " + bound.str();
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
    const auto write = [&](Spin component, const std::string& suffix,
                           const DensityProperty& property) {
        auto values =
            scalar_values(results.get(density_property(component, property.kind)), property.kind);
        const auto filename = prefix + "_" + suffix +
                              (property.kind == PropertyKind::Density ? "" : "_" + property.name) +
                              ".cube";
        const auto method = method_label(property.kind, results.backend_used(), esp_tolerance);
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
    if (o.ply)
        write_surface_ply(o, surface.mesh,
                          {{"esp", "Molecular electrostatic potential (hartree/e)"}},
                          {{"esp", surface.potential}});
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
    if (o.ply)
        write_surface_ply(o, surface.mesh, {{"alie", "Average local ionization energy (hartree)"}},
                          {{"alie", surface.ionization_energy}});
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

// Density components and properties selected with --spin and --property, and
// the contraction selected with --backend.
struct DensitySelection {
    std::vector<std::pair<Spin, std::string>> components; // with Cube filename suffix
    std::vector<DensityProperty> properties;
    chem::kumo::Backend backend = chem::kumo::Backend::AUTO;
};

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

std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (true) {
        const auto end = text.find(separator, start);
        parts.push_back(
            text.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) return parts;
        start = end + 1;
    }
}

// A location given as X,Y,Z in bohr or as @N, the position of atom N (1-based).
Eigen::Vector3d location(const std::string& spec, const Molecule& mol) {
    const auto& atoms = mol.atoms();
    if (!spec.empty() && spec[0] == '@') {
        const auto digits = spec.substr(1);
        if (digits.empty() || digits.size() > 9 ||
            !std::all_of(digits.begin(), digits.end(),
                         [](unsigned char c) { return std::isdigit(c); }))
            throw std::invalid_argument("invalid atom reference '" + spec + "' (use @N)");
        const auto n = std::stoul(digits);
        if (n == 0 || n > atoms.size())
            throw std::invalid_argument("atom " + spec + " does not exist (atoms are @1 to @" +
                                        std::to_string(atoms.size()) + ")");
        return atoms[n - 1].position;
    }
    const auto parts = split(spec, ',');
    Eigen::Vector3d point;
    try {
        if (parts.size() != 3) throw std::invalid_argument("three coordinates");
        for (int k = 0; k < 3; ++k)
            point[k] = real(parts[static_cast<std::size_t>(k)]);
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid location '" + spec +
                                    "' (use X,Y,Z in bohr or @N for atom N)");
    }
    return point;
}

// The locations of --line A:B or --plane A:B:C.
std::vector<Eigen::Vector3d> locations(const std::string& option, const std::string& spec,
                                       std::size_t count, const Molecule& mol) {
    const auto parts = split(spec, ':');
    if (parts.size() != count)
        throw std::invalid_argument(option + " needs " + std::to_string(count) +
                                    " locations separated by ':', for example " +
                                    (count == 2 ? "@1:@2" : "@1:@2:@3"));
    std::vector<Eigen::Vector3d> points;
    for (const auto& part : parts)
        points.push_back(location(part, mol));
    return points;
}

std::string coordinates(const Eigen::Vector3d& p) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << "(" << p.x() + 0.0 << ", " << p.y() + 0.0 << ", "
        << p.z() + 0.0 << ")";
    return out.str();
}

// A table of density properties: '#' comment lines describing the sampling
// and each column, a '#' line with the column names, then one row per point
// with the coordinates (bohr, fixed) and the values (scientific, 10 digits).
class Table {
public:
    Table(const DensitySelection& selection, std::vector<std::string> coordinate_columns)
        : selection_(selection), coordinate_columns_(std::move(coordinate_columns)) {
        for (const auto& [component, suffix] : selection.components)
            for (const auto& property : selection.properties)
                request_.add(density_property(component, property.kind));
    }

    const chem::kumo::Request& request() const noexcept { return request_; }

    void header(std::ostream& out, const std::vector<std::string>& description,
                chem::kumo::Backend used, double esp_tolerance) const {
        for (const auto& line : description)
            out << "# " << line << "\n";
        std::string names;
        for (const auto& name : coordinate_columns_)
            names += (names.empty() ? "" : " ") + name;
        for (const auto& [component, suffix] : selection_.components)
            for (const auto& property : selection_.properties) {
                const auto name = column(component, property);
                out << "# " << name << ": " << property.description << ", "
                    << component_label(component) << ", "
                    << method_label(property.kind, used, esp_tolerance) << "\n";
                names += " " + name;
            }
        out << "# " << names << "\n";
    }

    // The value columns of `results`, one vector per column.
    std::vector<std::vector<double>> columns(const chem::kumo::Results& results) const {
        std::vector<std::vector<double>> values;
        for (const auto& [component, suffix] : selection_.components)
            for (const auto& property : selection_.properties)
                values.push_back(scalar_values(
                    results.get(density_property(component, property.kind)), property.kind));
        return values;
    }

    // Row `q` of `columns`.
    static std::vector<double> at(const std::vector<std::vector<double>>& columns, std::size_t q) {
        std::vector<double> row;
        for (const auto& column : columns)
            row.push_back(column[q]);
        return row;
    }

    static void row(std::ostream& out, const std::vector<double>& coordinates,
                    const std::vector<double>& values) {
        std::ostringstream line;
        line << std::fixed << std::setprecision(8);
        // + 0.0 prints -0 as 0.
        for (std::size_t c = 0; c < coordinates.size(); ++c)
            line << (c > 0 ? " " : "") << std::setw(13) << coordinates[c] + 0.0;
        line << std::scientific << std::setprecision(10);
        for (const double value : values)
            line << " " << std::setw(17) << value;
        out << line.str() << "\n";
    }

private:
    static std::string column(Spin component, const DensityProperty& property) {
        const char* names[] = {"total", "alpha", "beta", "spin"};
        const auto index = component == Spin::Total   ? 0
                           : component == Spin::Alpha ? 1
                           : component == Spin::Beta  ? 2
                                                      : 3;
        return property.name + "_" + names[index];
    }

    const DensitySelection& selection_;
    std::vector<std::string> coordinate_columns_;
    chem::kumo::Request request_;
};

chem::kumo::EvaluationOptions evaluation_options(const DensitySelection& selection,
                                                 const Options& o, Console& console) {
    chem::kumo::EvaluationOptions options;
    options.backend = selection.backend;
    options.potential_tolerance = o.esp_tolerance;
    options.progress = console.progress();
    return options;
}

std::string sampling_summary(const DensitySelection& selection, std::size_t points) {
    return evaluation_summary(selection.components, selection.properties) + " at " +
           Console::grouped(points) + " point" + (points > 1 ? "s" : "");
}

// --point: one table on standard output, one row per point in the order given.
void print_points(const Molecule& mol, const Options& o, const DensitySelection& selection,
                  Console& console) {
    std::vector<Eigen::Vector3d> points;
    for (const auto& spec : o.points)
        points.push_back(location(spec, mol));
    const Table table(selection, {"x", "y", "z"});
    console.begin(sampling_summary(selection, points.size()));
    Eigen::Matrix<double, 3, Eigen::Dynamic> columns(3, static_cast<Eigen::Index>(points.size()));
    for (std::size_t p = 0; p < points.size(); ++p)
        columns.col(static_cast<Eigen::Index>(p)) = points[p];
    const auto results = chem::kumo::Evaluator(mol, evaluation_options(selection, o, console))
                             .evaluate(columns, table.request());
    console.end("Evaluated in");
    const auto values = table.columns(results);
    table.header(std::cout,
                 {"kumo " CHEM_KUMO_VERSION ", " + o.input.filename().string() + ": " +
                      std::to_string(points.size()) + " point" + (points.size() > 1 ? "s" : ""),
                  "x y z: position (bohr)"},
                 results.backend_used(), o.esp_tolerance);
    for (std::size_t p = 0; p < points.size(); ++p)
        Table::row(std::cout, {points[p].x(), points[p].y(), points[p].z()}, Table::at(values, p));
}

// --line and --plane: evaluate on the grid and write PREFIX_line.dat or
// PREFIX_plane.dat. Line rows give t, the signed distance from A along A->B;
// plane rows give u and v, the coordinates in the plane from A along the
// grid axes (u along A->B), with a blank line after each row of constant u,
// as gnuplot expects for surfaces.
void write_profile(const Molecule& mol, const Options& o, const DensitySelection& selection,
                   Console& console) {
    const bool is_line = !o.line.empty();
    const auto ends =
        is_line ? locations("--line", o.line, 2, mol) : locations("--plane", o.plane, 3, mol);
    const auto grid = is_line
                          ? chem::kumo::Grid::line(ends[0], ends[1], o.margin, o.spacing)
                          : chem::kumo::Grid::plane(ends[0], ends[1], ends[2], o.margin, o.spacing);
    const auto nu = static_cast<std::size_t>(grid.counts[0]);
    const auto nv = static_cast<std::size_t>(grid.counts[1]);
    const Eigen::Vector3d u = grid.steps.col(0).normalized(), v = grid.steps.col(1).normalized();

    std::vector<std::string> description;
    std::ostringstream geometry;
    geometry << std::setprecision(6);
    if (is_line)
        geometry << "Line from A " << coordinates(ends[0]) << " to B " << coordinates(ends[1])
                 << " bohr, extended by " << o.margin << " bohr: " << nu << " points "
                 << grid.steps.col(0).norm() << " bohr apart";
    else
        geometry << "Plane of A " << coordinates(ends[0]) << ", B " << coordinates(ends[1])
                 << " and C " << coordinates(ends[2]) << " bohr, margin " << o.margin
                 << " bohr: " << nu << " x " << nv << " points " << grid.steps.col(0).norm()
                 << " x " << grid.steps.col(1).norm() << " bohr apart";
    console.say(geometry.str());
    description.push_back("kumo " CHEM_KUMO_VERSION ", " + o.input.filename().string());
    description.push_back(geometry.str());
    if (is_line)
        description.push_back("t: distance from A towards B (bohr); x y z: position (bohr)");
    else {
        std::ostringstream axes;
        axes << std::setprecision(6) << "u, v: in-plane coordinates from A (bohr), u along A->B; "
             << "B at (" << (ends[1] - ends[0]).dot(u) << ", 0), C at ("
             << (ends[2] - ends[0]).dot(u) << ", " << (ends[2] - ends[0]).dot(v)
             << "); unit normal " << coordinates(grid.steps.col(2)) << "; x y z: position (bohr)";
        description.push_back(axes.str());
    }

    const Table table(selection, is_line ? std::vector<std::string>{"t", "x", "y", "z"}
                                         : std::vector<std::string>{"u", "v", "x", "y", "z"});
    console.begin(sampling_summary(selection, nu * nv));
    const auto results = chem::kumo::Evaluator(mol, evaluation_options(selection, o, console))
                             .evaluate(grid, table.request());
    console.end("Evaluated in");

    const auto filename = o.prefix + (is_line ? "_line.dat" : "_plane.dat");
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("cannot write " + filename);
    table.header(out, description, results.backend_used(), o.esp_tolerance);
    const auto values = table.columns(results);
    for (std::size_t i = 0; i < nu; ++i) {
        if (i > 0 && !is_line) out << "\n";
        for (std::size_t j = 0; j < nv; ++j) {
            const Eigen::Vector3d r = grid.origin + grid.steps.col(0) * static_cast<double>(i) +
                                      grid.steps.col(1) * static_cast<double>(j);
            const Eigen::Vector3d d = r - ends[0];
            std::vector<double> coordinates_row;
            if (is_line)
                coordinates_row = {d.dot(u), r.x(), r.y(), r.z()};
            else
                coordinates_row = {d.dot(u), d.dot(v), r.x(), r.y(), r.z()};
            Table::row(out, coordinates_row, Table::at(values, i * nv + j));
        }
    }
    out.flush();
    if (!out) throw std::runtime_error("write failed: " + filename);
    std::cout << filename << "\n";
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

// --surface-map: statistics of each selected property on the density
// isosurface, then (with --ply) the surface file. At most `listed` local
// extrema of each kind are printed.
void report_surface_map(const Molecule& mol, const Options& o, const DensitySelection& selection,
                        Console& console) {
    constexpr std::size_t listed = 20;
    std::vector<chem::kumo::Property> properties;
    std::vector<std::pair<std::string, std::string>> labels; // column name, description
    for (const auto& [component, suffix] : selection.components)
        for (const auto& property : selection.properties) {
            if (property.kind == PropertyKind::Gradient)
                throw std::invalid_argument(
                    "--surface-map cannot map the gradient (a vector); choose another property");
            properties.push_back(density_property(component, property.kind));
            const char* names[] = {"total", "alpha", "beta", "spin"};
            const auto index = component == Spin::Total   ? 0
                               : component == Spin::Alpha ? 1
                               : component == Spin::Beta  ? 2
                                                          : 3;
            auto name = property.name + "_" + names[index];
            std::replace(name.begin(), name.end(), '-', '_');
            labels.emplace_back(name, property.description + ", " + component_label(component));
        }
    chem::kumo::EvaluationOptions options;
    options.backend = selection.backend;
    options.potential_tolerance = o.esp_tolerance;
    console.begin(evaluation_summary(selection.components, selection.properties) +
                  " on the density isosurface");
    const auto surface =
        chem::kumo::map_on_density_isosurface(mol, properties, o.isovalue, o.spacing, options);
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
            << method_label(properties[f].kind, surface.backend_used, o.esp_tolerance) << "\n\n";
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
        if (o.surface_map) {
            const auto selection = select_density(o, mol, console);
            report_surface_map(mol, o, selection, console);
            console.finish();
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
        if (!o.points.empty() || !o.line.empty() || !o.plane.empty()) {
            const auto selection = select_density(o, mol, console);
            if (!o.points.empty())
                print_points(mol, o, selection, console);
            else
                write_profile(mol, o, selection, console);
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
            const auto selection = select_density(o, mol, console);
            write_density_properties(mol, grid, selection.components, selection.properties,
                                     selection.backend, o.esp_tolerance, o.prefix, console);
        }
        console.finish();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "kumo: " << e.what() << "\nTry 'kumo --help'.\n";
        return 1;
    }
}
