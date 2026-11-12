// Derived density properties against Multiwfn reference Cube files
// (provenance in tests/REFERENCE_CUBES.md). The Evaluator is run on the grid
// read from each Cube header and compared point by point.
//
// The Cube files store six significant digits, so rounding alone allows a
// relative difference of 5e-6. The tolerance, 2e-5 relative, also covers the
// largest differences observed: 1.2e-5 for the reduced gradient near the O-H
// bond critical point, where |grad rho| nearly cancels, and 9e-6 for
// sign(lambda2) rho at the oxygen nucleus.
#include "chem/kumo.hpp"
#include "chem/kumo/cube.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <string>

namespace {
using namespace chem::kumo;

constexpr double relative_tolerance = 2e-5;

void compare(const std::filesystem::path& fchk, const std::filesystem::path& cube_path,
             Property key, const std::string& label) {
    const auto molecule = Molecule::from_fchk(fchk);
    const auto reference = read_cube(cube_path);
    const auto& grid = reference.field.grid();
    const auto expected = reference.field.values();

    kumo_test::check(reference.field.atoms().size() == molecule.atoms().size(),
                     label + ": the Cube has the atoms of the FCHK");
    for (std::size_t a = 0; a < std::min(reference.field.atoms().size(), molecule.atoms().size());
         ++a)
        kumo_test::check(
            (reference.field.atoms()[a].position - molecule.atoms()[a].position).norm() < 1e-5,
            label + ": atom " + std::to_string(a + 1) + " is at the FCHK position");

    Request request;
    request.add(key);
    const auto actual = Evaluator(molecule).evaluate(grid, request).get(key).component();
    if (!kumo_test::check(actual.size() == expected.size(),
                          label + ": the Evaluator returns one value per Cube point"))
        return;

    double worst = 0.0;
    std::size_t worst_point = 0;
    bool same_sign = true;
    for (std::size_t p = 0; p < expected.size(); ++p) {
        const double error =
            std::abs(actual[p] - expected[p]) / std::max(std::abs(expected[p]), 1e-300);
        if (error > worst) {
            worst = error;
            worst_point = p;
        }
        same_sign = same_sign && (actual[p] < 0.0) == (expected[p] < 0.0);
    }
    std::cout << label << ": " << expected.size() << " points, max relative difference " << worst
              << " (point " << worst_point << ": " << actual[worst_point] << " versus "
              << expected[worst_point] << ")\n";
    kumo_test::check(same_sign, label + ": every value has the sign of the reference");
    kumo_test::check(worst <= relative_tolerance,
                     label + ": agrees with Multiwfn within " + kumo_test::num(relative_tolerance) +
                         " relative (max " + kumo_test::num(worst) + ")");
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " WATER.fchk O2.fchk MULTIWFN_DIRECTORY\n";
        return 2;
    }
    const std::filesystem::path water = argv[1], o2 = argv[2], directory = argv[3];
    compare(water, directory / "h2o_elf.cube.gz", density::total.elf(), "water ELF");
    compare(water, directory / "h2o_rdg.cube.gz", density::total.reduced_gradient(),
            "water reduced gradient");
    compare(water, directory / "h2o_sign_lambda2_rho.cube.gz",
            density::total.sign_lambda2_density(), "water sign(lambda2) rho");
    // Open shell: the total ELF uses the spin-polarized form.
    compare(o2, directory / "o2_elf.cube.gz", density::total.elf(), "O2 ELF");
    return kumo_test::exit_code();
}
