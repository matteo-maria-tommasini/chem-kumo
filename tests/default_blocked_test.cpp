#include <chem/kumo.hpp>
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    using namespace chem::kumo;
    return kumo_test::run([&] {
        for (int n = 1; n < argc; ++n) {
            const std::string name = argv[n];
            auto molecule = Molecule::from_fchk(argv[n]);
            auto grid = Grid::around(molecule.atoms(), 2.0, 0.65);
            double max_abs = 0.0;
            for (auto spin : {DensityComponent::Total, DensityComponent::Alpha,
                              DensityComponent::Beta, DensityComponent::Difference}) {
                const auto reference = molecule.density_pointwise(grid, spin);
                const auto current = molecule.density(grid, spin);
                const auto explicit_block = molecule.densities_blocked(grid, {spin}, 128).at(spin);
                const auto selected = molecule.densities(grid, {spin}).at(spin);
                auto a = reference.values(), b = current.values(), c = explicit_block.values(),
                     d = selected.values();
                const bool finite =
                    std::all_of(b.begin(), b.end(), [](double v) { return std::isfinite(v); }) &&
                    std::all_of(c.begin(), c.end(), [](double v) { return std::isfinite(v); }) &&
                    std::all_of(d.begin(), d.end(), [](double v) { return std::isfinite(v); });
                kumo_test::check(finite, name + " spin " + std::to_string(static_cast<int>(spin)) +
                                             ": blocked densities are finite");
                max_abs = std::max({max_abs, kumo_test::max_abs_diff(b, a, name + " density()"),
                                    kumo_test::max_abs_diff(c, a, name + " densities_blocked()"),
                                    kumo_test::max_abs_diff(d, a, name + " densities()")});
            }
            const auto combined = molecule.spin_densities(grid);
            for (auto spin : {DensityComponent::Total, DensityComponent::Alpha,
                              DensityComponent::Beta, DensityComponent::Difference}) {
                const auto reference_field = molecule.density_pointwise(grid, spin);
                const auto& field = spin == DensityComponent::Total   ? combined.total
                                    : spin == DensityComponent::Alpha ? combined.alpha
                                    : spin == DensityComponent::Beta  ? combined.beta
                                                                      : combined.spin;
                max_abs = std::max(max_abs,
                                   kumo_test::max_abs_diff(field.values(), reference_field.values(),
                                                           name + " spin_densities()"));
            }
            kumo_test::check(!(max_abs > 1e-9),
                             name +
                                 ": blocked densities match density_pointwise() within 1e-9 "
                                 "(max_abs=" +
                                 kumo_test::num(max_abs) + ")");
            std::cout << argv[n] << " max_abs=" << max_abs << '\n';
        }
    });
}
