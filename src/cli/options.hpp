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

/// Command-line options only: numerical evaluation stays in libkumo.
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

/// Print the help text to standard output.
void usage();

/// Parse and check the command line (without the program name handling of
/// --help and --version); throws std::invalid_argument with a message for
/// the user on any error.
Options parse(int argc, char** argv);

/// A real number, the whole string, finite; throws std::invalid_argument.
double real(const std::string& s);
} // namespace kumo_cli
