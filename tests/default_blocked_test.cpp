#include <chem/kumo.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    using namespace chem::kumo;
    for (int n = 1; n < argc; ++n) {
        auto molecule = Molecule::from_fchk(argv[n]);
        auto grid = Grid::around(molecule.atoms(), 2.0, 0.65);
        double max_abs = 0.0;
        for (auto spin : {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference}) {
            const auto reference = molecule.density_pointwise(grid, spin);
            const auto current = molecule.density(grid, spin);
            const auto explicit_block = molecule.densities_blocked(grid, {spin}, 128).at(spin);
            const auto selected = molecule.densities(grid, {spin}).at(spin);
            auto a = reference.values(), b = current.values(), c = explicit_block.values(), d = selected.values();
            for (std::size_t i = 0; i < a.size(); ++i) {
                if (!std::isfinite(b[i]) || !std::isfinite(c[i]) || !std::isfinite(d[i]))
                    throw std::runtime_error("nonfinite blocked density");
                max_abs = std::max({max_abs, std::abs(a[i]-b[i]),
                                    std::abs(a[i]-c[i]), std::abs(a[i]-d[i])});
            }
        }
        const auto combined = molecule.spin_densities(grid);
        for (auto spin : {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference}) {
            const auto reference_field = molecule.density_pointwise(grid, spin);
            const auto expected = reference_field.values();
            const auto& field = spin == Spin::Total ? combined.total :
                                spin == Spin::Alpha ? combined.alpha :
                                spin == Spin::Beta ? combined.beta : combined.spin;
            auto actual = field.values();
            for (std::size_t i = 0; i < actual.size(); ++i)
                max_abs = std::max(max_abs, std::abs(expected[i]-actual[i]));
        }
        if (max_abs > 1e-9) throw std::runtime_error("default blocked regression");
        std::cout << argv[n] << " max_abs=" << max_abs << '\n';
    }
}
