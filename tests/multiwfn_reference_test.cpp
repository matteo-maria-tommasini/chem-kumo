// Derived density properties and the electrostatic potential against
// Multiwfn reference Cube files (provenance in tests/REFERENCE_CUBES.md).
// The Evaluator is run on the grid read from each Cube header and compared
// point by point.
//
// The Cube files store six significant digits, so rounding alone allows a
// relative difference of 5e-6. The tolerance, 2e-5 relative, also covers the
// largest differences observed: 1.2e-5 for the reduced gradient near the O-H
// bond critical point, where |grad rho| nearly cancels, and 9e-6 for
// sign(lambda2) rho at the oxygen nucleus.
//
// The electrostatic potential also has an absolute tolerance of 1e-7
// hartree/e, for small values: the grid points are read from the Cube header,
// whose six decimals move them by up to 5e-7 bohr, and the potential changes
// by up to about 5e-8 over that distance away from the nuclei. Points within
// 1e-3 bohr of a nucleus are skipped: there the two programs treat the
// singular nuclear term differently (Multiwfn writes about 2e16 on the O2
// nucleus that its grid contains).
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

struct Tolerance {
    double relative = relative_tolerance;
    double absolute = 0.0;          ///< added to relative * |expected|
    double nucleus_exclusion = 0.0; ///< bohr; points this close to a nucleus are skipped
};

void compare(const std::filesystem::path& fchk, const std::filesystem::path& cube_path,
             Property key, const std::string& label, Tolerance tolerance = {}) {
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
    // Keep the Results alive: component() is a view into it.
    const Results results = Evaluator(molecule).evaluate(grid, request);
    const auto actual = results.get(key).component();
    if (!kumo_test::check(actual.size() == expected.size(),
                          label + ": the Evaluator returns one value per Cube point"))
        return;

    // Largest difference as a fraction of the allowed one.
    double worst = 0.0, worst_relative = 0.0;
    std::size_t worst_point = 0, skipped = 0;
    bool same_sign = true;
    const auto& atoms = molecule.atoms();
    for (std::size_t p = 0; p < expected.size(); ++p) {
        const std::size_t k = p % static_cast<std::size_t>(grid.counts[2]);
        const std::size_t j = (p / static_cast<std::size_t>(grid.counts[2])) %
                              static_cast<std::size_t>(grid.counts[1]);
        const std::size_t i =
            p / static_cast<std::size_t>(grid.counts[2]) / static_cast<std::size_t>(grid.counts[1]);
        const Eigen::Vector3d point =
            grid.origin + grid.steps * Eigen::Vector3d(double(i), double(j), double(k));
        if (std::any_of(atoms.begin(), atoms.end(), [&](const Atom& atom) {
                return (atom.position - point).norm() < tolerance.nucleus_exclusion;
            })) {
            ++skipped;
            continue;
        }
        const double difference = std::abs(actual[p] - expected[p]);
        const double allowed = tolerance.relative * std::abs(expected[p]) + tolerance.absolute;
        const double ratio = difference / std::max(allowed, 1e-300);
        if (ratio > worst) {
            worst = ratio;
            worst_point = p;
        }
        if (std::abs(expected[p]) > 1e3 * tolerance.absolute)
            worst_relative = std::max(worst_relative, difference / std::abs(expected[p]));
        if (std::abs(expected[p]) > tolerance.absolute)
            same_sign = same_sign && (actual[p] < 0.0) == (expected[p] < 0.0);
    }
    std::cout << label << ": " << expected.size() - skipped << " points (" << skipped
              << " skipped near nuclei), max difference " << worst << " of the tolerance (point "
              << worst_point << ": " << actual[worst_point] << " versus " << expected[worst_point]
              << "), max relative difference " << worst_relative << '\n';
    kumo_test::check(same_sign, label + ": every value has the sign of the reference");
    kumo_test::check(worst <= 1.0, label + ": agrees with Multiwfn within " +
                                       kumo_test::num(tolerance.relative) + " relative + " +
                                       kumo_test::num(tolerance.absolute) + " absolute (max " +
                                       kumo_test::num(worst) + " of the tolerance)");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 3, 3, "WATER.fchk O2.fchk MULTIWFN_DIRECTORY"))
        return 2;
    return kumo_test::run([&] {
        const std::filesystem::path water = argv[1], o2 = argv[2], directory = argv[3];
        compare(water, directory / "h2o_elf.cube.gz", density::total.elf(), "water ELF");
        compare(water, directory / "h2o_rdg.cube.gz", density::total.reduced_gradient(),
                "water reduced gradient");
        compare(water, directory / "h2o_sign_lambda2_rho.cube.gz",
                density::total.sign_lambda2_density(), "water sign(lambda2) rho");
        // Open shell: the total ELF uses the spin-polarized form.
        compare(o2, directory / "o2_elf.cube.gz", density::total.elf(), "O2 ELF");
        const Tolerance potential{relative_tolerance, 1e-7, 1e-3};
        compare(water, directory / "h2o_esp.cube.gz", density::total.electrostatic_potential(),
                "water electrostatic potential", potential);
        compare(o2, directory / "o2_esp.cube.gz", density::total.electrostatic_potential(),
                "O2 electrostatic potential", potential);
    });
}
