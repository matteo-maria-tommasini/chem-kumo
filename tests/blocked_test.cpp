#include "chem/kumo.hpp"
#include "molecule_test_access.hpp"
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
            auto mol = Molecule::from_fchk(argv[n]);
            auto grid = Grid::around(mol.atoms(), 2.0, 0.65);
            double maxerr = 0;
            for (auto size : {std::size_t(1), std::size_t(16), std::size_t(32), std::size_t(64),
                              std::size_t(128), std::size_t(256)}) {
                auto ref = kumo_test::Kernels::densities_blocked(
                    mol, grid,
                    {DensityComponent::Total, DensityComponent::Alpha, DensityComponent::Beta,
                     DensityComponent::Difference});
                auto blk = kumo_test::Kernels::densities_blocked(
                    mol, grid,
                    {DensityComponent::Total, DensityComponent::Alpha, DensityComponent::Beta,
                     DensityComponent::Difference},
                    size);
                for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                                       DensityComponent::Beta, DensityComponent::Difference})
                    maxerr =
                        std::max(maxerr, kumo_test::max_abs_diff(
                                             blk.at(component).values(), ref.at(component).values(),
                                             name + " block " + std::to_string(size)));
            }
            auto only = kumo_test::Kernels::densities_blocked(
                mol, grid, {DensityComponent::Total, DensityComponent::Total}, 32);
            kumo_test::check(only.size() == 1,
                             name +
                                 ": a duplicate component is evaluated once (expected 1 field, "
                                 "got " +
                                 std::to_string(only.size()) + ")");
            kumo_test::check_throws<std::invalid_argument>(
                [&] {
                    (void)kumo_test::Kernels::densities_blocked(mol, grid,
                                                                {DensityComponent::Total}, 0);
                },
                name + ": densities_blocked rejects block size 0");
            kumo_test::check(!(maxerr > 1e-9),
                             name + ": blocked densities match densities() within 1e-9 (max_abs=" +
                                 kumo_test::num(maxerr) + ")");
            std::cout << argv[n] << " max_abs=" << maxerr << "\n";
        }
    });
}
