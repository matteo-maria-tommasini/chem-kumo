#include "chem/kumo.hpp"
#include <Eigen/Core>
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
        Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 137);
        for (Eigen::Index j=0; j<points.cols(); ++j)
            points.col(j) = Eigen::Vector3d(0.21+0.043*j, -0.37+0.017*j, 0.14-0.031*j);
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(1);
#endif
        for (const auto spin : {chem::kumo::Spin::Alpha, chem::kumo::Spin::Beta}) {
            const auto n = mol.orbital_count(spin);
            if (n == 0) throw std::runtime_error("missing molecular orbitals");
            const std::vector<std::size_t> indices = {1, std::min<std::size_t>(3,n),
                std::min<std::size_t>(static_cast<std::size_t>(mol.occupied_count(spin))+1,n)};
            const auto reference = mol.orbitals(grid, indices, spin);
            const auto derivative_reference = mol.orbital_derivatives_blocked(points, indices[1], spin, 64);
            for (int threads : {1, 2, 4}) {
#ifdef _OPENMP
                omp_set_num_threads(threads);
#else
                (void)threads;
#endif
                const auto fields = mol.orbitals(grid, indices, spin);
                if (fields.size() != reference.size()) throw std::runtime_error("MO count mismatch");
                for (std::size_t k=0; k<fields.size(); ++k) {
                    const auto a=reference[k].values(), b=fields[k].values();
                    if (a.size()!=b.size()) throw std::runtime_error("MO grid size mismatch");
                    for (std::size_t i=0; i<a.size(); ++i)
                        if (!std::isfinite(b[i]) || std::abs(a[i]-b[i]) > 1e-10)
                            throw std::runtime_error("MO parallel grid mismatch");
                }
                for (std::size_t block : {std::size_t(1), std::size_t(32), std::size_t(64), std::size_t(128)}) {
                    const auto derivatives = mol.orbital_derivatives_blocked(points, indices[1], spin, block);
                    if (derivatives.size()!=derivative_reference.size()) throw std::runtime_error("MO derivative size mismatch");
                    for (std::size_t j=0; j<derivatives.size(); ++j) {
                        const auto& a=derivative_reference[j];
                        const auto& b=derivatives[j];
                        if (!std::isfinite(b.value) || !b.gradient.allFinite() || !b.hessian.allFinite() ||
                            std::abs(a.value-b.value)>1e-10 ||
                            (a.gradient-b.gradient).cwiseAbs().maxCoeff()>1e-10 ||
                            (a.hessian-b.hessian).cwiseAbs().maxCoeff()>1e-10)
                            throw std::runtime_error("MO parallel derivative mismatch");
                    }
                }
            }
        }
        bool rejected=false;
        try { (void)mol.orbital_derivatives_blocked(points, 1, chem::kumo::Spin::Alpha, 0); }
        catch (const std::invalid_argument&) { rejected=true; }
        if (!rejected) throw std::runtime_error("zero block size accepted");
        std::cout << argv[file] << " parallel orbital grid/derivatives OK\n";
    }
}
