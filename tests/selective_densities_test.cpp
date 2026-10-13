#include "chem/kumo.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    using namespace chem::kumo;
    for (int n = 1; n < argc; ++n) {
        auto molecule = Molecule::from_fchk(argv[n]);
        auto grid = Grid::around(molecule.atoms(), 1.5, 0.6);
        auto all = molecule.spin_densities(grid);
        auto selected = molecule.densities(grid, {Spin::Total, Spin::Difference, Spin::Total});
        if (selected.size() != 2 || !molecule.densities(grid, {}).empty()) return 1;
        double error = 0;
        for (auto [key, expected] : {std::pair{Spin::Total, &all.total},
                                      std::pair{Spin::Difference, &all.spin}}) {
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
        if (!std::isfinite(error) || error > 1e-9) return 1;
    }
    return 0;
}
