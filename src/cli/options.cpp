#include "cli/options.hpp"
#include <cctype>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace kumo_cli {
namespace {
std::string lower(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
} // namespace

double real(const std::string& s) {
    std::size_t n = 0;
    const auto v = std::stod(s, &n);
    if (n != s.size() || !std::isfinite(v))
        throw std::invalid_argument("invalid real number: " + s);
    return v;
}

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
    o.mode = o.list              ? Mode::ListOrbitals
             : o.surface_map     ? Mode::SurfaceMap
             : o.esp_surface     ? Mode::EspSurface
             : o.alie_surface    ? Mode::AlieSurface
             : !o.points.empty() ? Mode::Points
             : !o.line.empty()   ? Mode::Line
             : !o.plane.empty()  ? Mode::Plane
             : o.mo_mode         ? Mode::OrbitalCubes
                                 : Mode::DensityCubes;
    return o;
}
} // namespace kumo_cli
