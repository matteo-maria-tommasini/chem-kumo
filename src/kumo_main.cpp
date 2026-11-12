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
              << "  --property density,gradient,laplacian,tau,elf,rdg,sign-lambda2-rho\n"
              << "                             Density properties to write (default density;\n"
              << "                             gradient is the norm |grad rho|; elf for total,\n"
              << "                             alpha or beta; rdg and sign-lambda2-rho for the\n"
              << "                             total density only)\n"
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
        } else if (a == "--spacing")
            o.spacing = real(value());
        else if (a == "--margin")
            o.margin = real(value());
        else if (a == "--output") {
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
};

// Parse --property: a comma-separated list, duplicates ignored, user order kept.
std::vector<DensityProperty> select_properties(const std::string& list) {
    const DensityProperty known[] = {
        {PropertyKind::Density, "density", "Electron density (e/bohr^3)"},
        {PropertyKind::Gradient, "gradient",
         "Electron density gradient norm |grad rho| (e/bohr^4)"},
        {PropertyKind::Laplacian, "laplacian", "Electron density Laplacian (e/bohr^5)"},
        {PropertyKind::Tau, "tau", "Positive kinetic energy density tau (hartree/bohr^3)"},
        {PropertyKind::Elf, "elf", "Electron localization function ELF (dimensionless)"},
        {PropertyKind::ReducedGradient, "rdg",
         "Reduced density gradient s (dimensionless, capped at 100)"},
        {PropertyKind::SignLambda2Density, "sign-lambda2-rho",
         "sign(lambda2) * electron density (e/bohr^3)"}};
    std::vector<DensityProperty> result;
    std::size_t start = 0;
    while (true) {
        const auto comma = list.find(',', start);
        const auto token =
            list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        const auto* match = std::find_if(std::begin(known), std::end(known),
                                         [&](const DensityProperty& p) { return p.name == token; });
        if (match == std::end(known))
            throw std::invalid_argument(
                "unknown --property '" + token +
                "' (use density, gradient, laplacian, tau, elf, rdg or sign-lambda2-rho)");
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
    default: return d.value();
    }
}

// ELF is not defined for the spin difference, and the reduced gradient and
// sign(lambda2) rho only for the total density.
void check_components(const std::vector<std::pair<Spin, std::string>>& components,
                      const std::vector<DensityProperty>& properties) {
    for (const auto& property : properties)
        for (const auto& component : components) {
            const Spin spin = component.first;
            if (property.kind == PropertyKind::Elf && spin == Spin::Difference)
                throw std::invalid_argument(
                    "--property elf is not defined for --spin difference (or all); "
                    "use total, alpha or beta");
            if ((property.kind == PropertyKind::ReducedGradient ||
                 property.kind == PropertyKind::SignLambda2Density) &&
                spin != Spin::Total)
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
        what.push_back(property.kind == PropertyKind::Density     ? "value"
                       : property.kind == PropertyKind::Gradient  ? "gradient norm"
                       : property.kind == PropertyKind::Laplacian ? "Laplacian"
                       : property.kind == PropertyKind::Tau       ? "kinetic energy density tau"
                       : property.kind == PropertyKind::Elf       ? "ELF"
                       : property.kind == PropertyKind::ReducedGradient ? "reduced density gradient"
                                                                        : "sign(lambda2) rho");
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
                              chem::kumo::Backend backend, const std::string& prefix,
                              Console& console) {
    chem::kumo::Request request;
    for (const auto& [component, suffix] : components)
        for (const auto& property : properties)
            request.add(density_property(component, property.kind));
    console.begin(evaluation_summary(components, properties));
    chem::kumo::EvaluationOptions options;
    options.backend = backend;
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
        if (o.list) {
            if (o.spin.empty() || o.spin == "alpha" || o.spin == "both") list_mos(mol, Spin::Alpha);
            if (o.spin == "beta" || o.spin == "both") list_mos(mol, Spin::Beta);
            if (!(o.spin.empty() || o.spin == "alpha" || o.spin == "beta" || o.spin == "both"))
                throw std::invalid_argument("--list-mo accepts alpha, beta, or both");
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
            const auto backend = select_backend(o.backend, mol);
            console.say(backend_summary(backend, mol));
            write_density_properties(mol, grid, components, properties, backend, o.prefix, console);
        }
        console.finish();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "kumo: " << e.what() << "\nTry 'kumo --help'.\n";
        return 1;
    }
}
