#pragma once
// Options of the kumo command-line program.
#include <filesystem>
#include <string>
#include <vector>

namespace kumo_cli {
namespace fs = std::filesystem;

/// What the program does, decided by the options (one mode per run).
enum class Mode {
    DensityCubes, ///< default: Cube files of density properties
    OrbitalCubes, ///< --mo
    ListOrbitals, ///< --list-mo
    EspSurface,   ///< --esp-surface
    AlieSurface,  ///< --alie-surface
    SurfaceMap,   ///< --surface-map
    Points,       ///< --point
    Line,         ///< --line
    Plane,        ///< --plane
};

/// The parsed command line; all computation is left to libkumo. String
/// options hold the option's argument (lower-cased for --mo, --spin,
/// --property and --backend) and are empty when the option is not given.
/// The *_set flags record which numerical options were given, because the
/// defaults depend on the mode; parse() fills in those defaults.
struct Options {
    Mode mode = Mode::DensityCubes;
    fs::path input;
    std::string selector;
    std::string spin;
    std::string properties;
    std::string backend;
    std::string prefix;
    double spacing = 0.20;
    double margin = 5.0;
    double esp_tolerance = 1e-8; // EvaluationOptions::potential_tolerance default
    bool esp_tolerance_set = false;
    double isovalue = 1e-3; // density isosurface of the surface modes, e/bohr^3
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

/// Print the help text to standard output.
void usage();

/// Parse and check the command-line arguments argv[1] to argv[argc - 1],
/// and set the mode. main() handles a lone --help or --version before
/// calling this. Throws std::invalid_argument with a message for the user on
/// any error (std::out_of_range for a number too large for a double).
Options parse(int argc, char** argv);

/// Parse `s`, which must be entirely a finite real number. Throws
/// std::invalid_argument otherwise (std::out_of_range, from std::stod, for a
/// value outside the range of double).
double real(const std::string& s);
} // namespace kumo_cli
