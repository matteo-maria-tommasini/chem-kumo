// densities() with a subset of components (duplicates allowed, empty list
// allowed) must match spin_densities() and single-component density().
#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    using namespace chem::kumo;
    return kumo_test::run([&] {
        for (int n = 1; n < argc; ++n) {
            auto molecule = Molecule::from_fchk(argv[n]);
            auto grid = Grid::around(molecule.atoms(), 1.5, 0.6);
            auto all = molecule.spin_densities(grid);
            auto selected = molecule.densities(grid, {Spin::Total, Spin::Difference, Spin::Total});
            kumo_test::check(selected.size() == 2,
                             "duplicate components are merged (expected 2 fields, got " +
                                 std::to_string(selected.size()) + ")");
            kumo_test::check(molecule.densities(grid, {}).empty(),
                             "an empty component list returns no fields");
            double error = 0;
            for (auto [key, expected] :
                 {std::pair{Spin::Total, &all.total}, std::pair{Spin::Difference, &all.spin}}) {
                error = std::max(error, kumo_test::max_abs_diff(selected.at(key).values(),
                                                                expected->values(),
                                                                "densities() vs spin_densities()"));
            }
            for (Spin key : {Spin::Alpha, Spin::Beta, Spin::Total, Spin::Difference}) {
                auto one = molecule.densities(grid, {key});
                auto field = molecule.density(grid, key);
                error =
                    std::max(error, kumo_test::max_abs_diff(one.at(key).values(), field.values(),
                                                            "densities() vs density()"));
            }
            std::cout << argv[n] << " max_abs=" << error << '\n';
            kumo_test::check(std::isfinite(error) && error <= 1e-9,
                             std::string(argv[n]) +
                                 ": selected densities match within 1e-9 (max_abs=" +
                                 kumo_test::num(error) + ")");
        }
    });
}
