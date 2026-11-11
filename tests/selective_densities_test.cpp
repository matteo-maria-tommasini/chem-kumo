// densities() with a subset of components (duplicates allowed, empty list
// allowed) must match spin_densities() and single-component density().
#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " FILE1.fchk FILE2.fchk\n";
        return 2;
    }
    using namespace chem::kumo;
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
            const auto actual = selected.at(key).values();
            for (std::size_t i = 0; i < actual.size(); ++i)
                error = std::max(error, std::abs(actual[i] - expected->values()[i]));
        }
        for (Spin key : {Spin::Alpha, Spin::Beta, Spin::Total, Spin::Difference}) {
            auto one = molecule.densities(grid, {key});
            auto field = molecule.density(grid, key);
            for (std::size_t i = 0; i < field.values().size(); ++i)
                error = std::max(error, std::abs(one.at(key).values()[i] - field.values()[i]));
        }
        std::cout << argv[n] << " max_abs=" << error << '\n';
        kumo_test::check(std::isfinite(error) && error <= 1e-9,
                         std::string(argv[n]) + ": selected densities match within 1e-9 (max_abs=" +
                             kumo_test::num(error) + ")");
    }
    return kumo_test::exit_code();
}
