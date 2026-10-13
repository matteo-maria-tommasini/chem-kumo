#include "chem/kumo.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    try {
        for (int t = 1; t < argc; ++t) {
            const auto molecule = chem::kumo::Molecule::from_fchk(argv[t]);
            auto grid = chem::kumo::Grid::around(molecule.atoms(), 1.5, 0.6);
            auto all = molecule.spin_densities(grid);
            const auto a = molecule.density(grid, chem::kumo::Spin::Alpha);
            const auto b = molecule.density(grid, chem::kumo::Spin::Beta);
            const auto total = molecule.density(grid, chem::kumo::Spin::Total);
            const auto spin = molecule.density(grid, chem::kumo::Spin::Difference);
            double worst = 0.0;
            for (std::size_t i = 0; i < all.total.values().size(); ++i) {
                const double av = all.alpha.values()[i], bv = all.beta.values()[i];
                worst = std::max({worst, std::abs(av-a.values()[i]),
                                  std::abs(bv-b.values()[i]),
                                  std::abs(all.total.values()[i]-total.values()[i]),
                                  std::abs(all.spin.values()[i]-spin.values()[i]),
                                  std::abs(all.total.values()[i]-(av+bv)),
                                  std::abs(all.spin.values()[i]-(av-bv))});
            }
            std::cout << argv[t] << " points=" << all.total.values().size()
                      << " max_abs=" << worst << '\n';
            if (!std::isfinite(worst) || worst > 1e-9) return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
