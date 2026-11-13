// Electrostatic potential surface analysis against Multiwfn.
//
// The references (tests/data/multiwfn/*_esp_surface.txt, provenance in
// tests/REFERENCE_CUBES.md) are Multiwfn's quantitative molecular surface
// analysis of the ESP on the rho = 0.001 isosurface with a grid spacing of
// 0.05 bohr, near its converged limit. kumo runs with its default spacing of
// 0.25 bohr, so the differences are kumo's discretization error at that
// spacing, which falls as the spacing squared, plus Multiwfn's at 0.05 bohr.
// The largest observed differences are 1.6e-3 relative for areas and
// averages (the small potential of O2), 3.4e-3 for variances (O2, where kumo
// itself moves by 5e-3 from 0.25 to 0.05 bohr), 1.4e-5 hartree/e for extremum
// values and 0.13 bohr for extremum positions; the tolerances allow about
// three times as much.

#include "chem/kumo.hpp"
#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Multiwfn's conversion factors (define.f90): its kcal/mol values are
// converted back to hartree with the same factor.
constexpr double multiwfn_kcal = 627.51;
constexpr double bohr_in_angstrom = 0.529177210903;

std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);)
        lines.push_back(line);
    return lines;
}

// Numbers among the whitespace-separated tokens of `text`, parentheses and
// asterisks removed.
std::vector<double> numbers(std::string text) {
    std::replace_if(
        text.begin(), text.end(), [](char c) { return c == '(' || c == ')' || c == '*'; }, ' ');
    std::istringstream in(text);
    std::vector<double> result;
    for (std::string token; in >> token;) {
        try {
            std::size_t used = 0;
            const double value = std::stod(token, &used);
            if (used == token.size()) result.push_back(value);
        } catch (const std::exception&) {
        }
    }
    return result;
}

// The index-th number after `label` on the first line containing it.
double value_after(const std::vector<std::string>& lines, const std::string& label,
                   std::size_t index) {
    for (const auto& line : lines) {
        const auto at = line.find(label);
        if (at == std::string::npos) continue;
        const auto values = numbers(line.substr(at + label.size()));
        if (index >= values.size())
            throw std::runtime_error("reference line '" + line + "' lacks a value");
        return values[index];
    }
    throw std::runtime_error("reference lacks '" + label + "'");
}

struct Extremum {
    double value = 0.0;       // hartree/e
    Eigen::Vector3d position; // bohr
};

// The table after "The number of surface minima:" (or maxima).
std::vector<Extremum> extrema(const std::vector<std::string>& lines, const std::string& kind) {
    const std::string label = "The number of surface " + kind + ":";
    for (std::size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find(label) == std::string::npos) continue;
        const auto count = static_cast<std::size_t>(value_after(lines, label, 0));
        std::vector<Extremum> result;
        for (std::size_t row = i + 2; row < i + 2 + count; ++row) {
            const auto v = numbers(lines.at(row)); // #, a.u., eV, kcal/mol, x, y, z (Angstrom)
            if (v.size() != 7) throw std::runtime_error("malformed extremum: " + lines.at(row));
            result.push_back({v[1], Eigen::Vector3d(v[4], v[5], v[6]) / bohr_in_angstrom});
        }
        return result;
    }
    throw std::runtime_error("reference lacks '" + label + "'");
}

// Largest ratio of |difference| to tolerance, for the report.
double worst_ratio = 0.0;

void compare(double actual, double expected, double rel, double abs, const std::string& label) {
    const double tolerance = std::max(abs, rel * std::abs(expected));
    worst_ratio = std::max(worst_ratio, std::abs(actual - expected) / tolerance);
    kumo_test::check_close(actual, expected, rel, abs, label);
}

// Each reference extremum has a kumo extremum of the same kind within
// `distance` bohr and `value` hartree/e, and the counts agree.
void compare_extrema(const chem::kumo::PotentialSurface& surface,
                     const std::vector<std::size_t>& found, const std::vector<Extremum>& expected,
                     const std::string& label) {
    if (!kumo_test::check(found.size() == expected.size(),
                          label + ": " + std::to_string(expected.size()) + " local extrema (got " +
                              std::to_string(found.size()) + ")"))
        return;
    std::vector<bool> used(found.size(), false);
    for (const auto& e : expected) {
        std::size_t best = found.size();
        double nearest = 0.0;
        for (std::size_t f = 0; f < found.size(); ++f) {
            const double d = (surface.mesh.vertices[found[f]] - e.position).norm();
            if (!used[f] && (best == found.size() || d < nearest)) {
                best = f;
                nearest = d;
            }
        }
        used[best] = true;
        const std::string what = label + " at (" + kumo_test::num(e.position.x()) + ", " +
                                 kumo_test::num(e.position.y()) + ", " +
                                 kumo_test::num(e.position.z()) + ") bohr";
        // A surface vertex lies within about a grid spacing of the extremum.
        compare(nearest, 0.0, 0.0, 0.4, what + ": distance to the nearest kumo extremum");
        compare(surface.potential[found[best]], e.value, 0.0, 4e-5, what + ": value");
    }
}

void check_case(const std::string& fchk, const std::string& reference, bool local_extrema) {
    const auto lines = read_lines(reference);
    const auto mol = chem::kumo::Molecule::from_fchk(fchk);
    const auto surface = chem::kumo::electrostatic_potential_surface(mol);
    const auto& s = surface.statistics;
    const std::string label = reference;
    const auto kcal = [&](const std::string& key, std::size_t index) {
        return value_after(lines, key, index) / multiwfn_kcal;
    };
    const auto kcal2 = [&](const std::string& key) {
        return value_after(lines, key, 1) / (multiwfn_kcal * multiwfn_kcal);
    };

    const double area_rel = 5e-3;
    compare(surface.mesh.volume(), value_after(lines, "Volume:", 0), area_rel, 0,
            label + ": volume");
    compare(s.area, value_after(lines, "Overall surface area:", 0), area_rel, 0, label + ": area");
    compare(s.positive_area, value_after(lines, "Positive surface area:", 0), area_rel, 0,
            label + ": positive area");
    compare(s.negative_area, value_after(lines, "Negative surface area:", 0), area_rel, 0,
            label + ": negative area");
    compare(s.nonpolar_area / s.area,
            value_after(lines, "Nonpolar surface area (|ESP| <= 10 kcal/mol):", 1) / 100, 0, 2e-3,
            label + ": nonpolar area fraction");

    compare(s.minimum, kcal("Minimal value:", 0), 0, 4e-5, label + ": global minimum");
    compare(s.maximum, kcal("Minimal value:", 1), 0, 4e-5, label + ": global maximum");

    // Averages and deviations: relative tolerances, with a floor of 2e-6
    // hartree/e for the nearly neutral surface of O2.
    const double rel = 5e-3, abs = 2e-6;
    compare(s.average, value_after(lines, "Overall average value:", 0), rel, abs,
            label + ": average");
    compare(s.positive_average, value_after(lines, "Positive average value:", 0), rel, abs,
            label + ": positive average");
    compare(s.negative_average, value_after(lines, "Negative average value:", 0), rel, abs,
            label + ": negative average");
    compare(s.charge_separation, kcal("Internal charge separation (Pi):", 1), rel, abs,
            label + ": Pi");
    compare(s.polarity_index, kcal("Molecular polarity index (MPI):", 1), rel, abs,
            label + ": MPI");
    const double variance_rel = 1e-2, variance_abs = 2e-8;
    compare(s.total_variance, kcal2("Overall variance (sigma^2_tot):"), variance_rel, variance_abs,
            label + ": total variance");
    compare(s.positive_variance, kcal2("Positive variance:"), variance_rel, variance_abs,
            label + ": positive variance");
    compare(s.negative_variance, kcal2("Negative variance:"), variance_rel, variance_abs,
            label + ": negative variance");
    compare(s.balance, value_after(lines, "Balance of charges (nu):", 0), 0, 1e-3,
            label + ": balance nu");

    if (local_extrema) {
        compare_extrema(surface, s.local_minima, extrema(lines, "minima"), label + " minimum");
        compare_extrema(surface, s.local_maxima, extrema(lines, "maxima"), label + " maximum");
    }
}

} // namespace

int main(int argc, char** argv) {
    // Arguments: FCHK REFERENCE yes|no (compare local extrema), repeated.
    if (!kumo_test::require_args(argc, argv, 3, kumo_test::unlimited,
                                 "FILE.fchk REFERENCE.txt yes|no ...") ||
        (argc - 1) % 3 != 0)
        return 2;
    return kumo_test::run([&] {
        for (int i = 1; i + 2 < argc; i += 3) {
            const std::string extrema_flag = argv[i + 2];
            if (extrema_flag != "yes" && extrema_flag != "no")
                throw std::invalid_argument("expected yes or no, got " + extrema_flag);
            check_case(argv[i], argv[i + 1], extrema_flag == "yes");
        }
        std::cout << "largest difference: " << kumo_test::num(worst_ratio) << " of its tolerance\n";
        if (kumo_test::failures == 0) std::cout << "ESP surface analysis agrees with Multiwfn\n";
    });
}
