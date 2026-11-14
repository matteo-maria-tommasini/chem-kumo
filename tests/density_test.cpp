// Density components on a coarse grid: the total density is finite and
// non-negative, total = alpha + beta and spin = alpha - beta at every point.
#include <chem/kumo.hpp>
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <cmath>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    try {
        for (int f = 1; f < 3; ++f) {
            const std::string file = argv[f];
            auto mol = chem::kumo::Molecule::from_fchk(file);
            auto g = chem::kumo::Grid::around(mol.atoms(), 2.0, 0.7);
            auto t = kumo_test::Kernels::density(mol, g, chem::kumo::DensityComponent::Total);
            auto a = kumo_test::Kernels::density(mol, g, chem::kumo::DensityComponent::Alpha);
            auto b = kumo_test::Kernels::density(mol, g, chem::kumo::DensityComponent::Beta);
            auto s = kumo_test::Kernels::density(mol, g, chem::kumo::DensityComponent::Difference);
            // Report only the first failing point of each file.
            for (std::size_t i = 0; i < t.values().size(); ++i) {
                const double total = t.values()[i], alpha = a.values()[i], beta = b.values()[i];
                const double spin = s.values()[i];
                const std::string at = file + " point " + std::to_string(i) + ": ";
                if (!kumo_test::check(std::isfinite(total) && total >= -1e-8,
                                      at + "total density finite and >= -1e-8 (got " +
                                          kumo_test::num(total) + ")") ||
                    !kumo_test::check(std::abs(total - alpha - beta) <= 1e-9,
                                      at + "total = alpha + beta within 1e-9 (difference " +
                                          kumo_test::num(total - alpha - beta) + ")") ||
                    !kumo_test::check(std::abs(spin - alpha + beta) <= 1e-9,
                                      at + "spin = alpha - beta within 1e-9 (difference " +
                                          kumo_test::num(spin - alpha + beta) + ")"))
                    break;
            }
            std::cout << file << " points=" << t.values().size() << " basis=" << mol.basis_size()
                      << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
    return kumo_test::exit_code();
}
