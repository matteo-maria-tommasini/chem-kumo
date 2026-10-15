#include "chem/kumo.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    for (int file = 1; file < argc; ++file) {
        const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
        const auto grid = chem::kumo::Grid::around(mol.atoms(), 2.0, 0.55);
        const auto components = {chem::kumo::Spin::Total, chem::kumo::Spin::Alpha,
                                 chem::kumo::Spin::Beta, chem::kumo::Spin::Difference};
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(1);
#endif
        const auto reference = mol.densities_blocked(grid, components, 64);
        for (int threads : {1, 2, 4}) {
#ifdef _OPENMP
            omp_set_num_threads(threads);
#else
            (void)threads;
#endif
            for (std::size_t block : {std::size_t(1), std::size_t(32), std::size_t(64), std::size_t(128)}) {
                const auto fields = mol.densities_blocked(grid, components, block);
                for (auto spin : components) {
                    const auto a = reference.at(spin).values();
                    const auto b = fields.at(spin).values();
                    if (a.size() != b.size()) throw std::runtime_error("parallel density size mismatch");
                    for (std::size_t i = 0; i < a.size(); ++i) {
                        if (!std::isfinite(b[i]) || std::abs(a[i] - b[i]) > 1e-10)
                            throw std::runtime_error("parallel density numerical mismatch");
                    }
                }
            }
        }
        std::cout << argv[file] << " parallel/serial density comparison OK\n";
    }
}
