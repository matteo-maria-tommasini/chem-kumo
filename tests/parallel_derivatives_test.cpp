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
        Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 137);
        for (Eigen::Index i = 0; i < points.cols(); ++i) {
            points.col(i) = Eigen::Vector3d(0.11 + 0.073*i, -1.3 + 0.041*i, 0.27 - 0.037*i);
        }
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(1);
#endif
        for (auto spin : {chem::kumo::Spin::Total, chem::kumo::Spin::Alpha,
                          chem::kumo::Spin::Beta, chem::kumo::Spin::Difference}) {
            const auto reference = mol.density_derivatives_blocked(points, spin, 32);
            for (int threads : {1, 2, 4}) {
#ifdef _OPENMP
                omp_set_num_threads(threads);
#else
                (void)threads;
#endif
                for (std::size_t block : {std::size_t(1), std::size_t(16), std::size_t(32), std::size_t(64)}) {
                    const auto actual = mol.density_derivatives_blocked(points, spin, block);
                    if (actual.size() != reference.size()) throw std::runtime_error("derivative size mismatch");
                    for (std::size_t i = 0; i < actual.size(); ++i) {
                        const auto& a = actual[i];
                        const auto& b = reference[i];
                        if (!std::isfinite(a.value) || !a.gradient.allFinite() || !a.hessian.allFinite() ||
                            !std::isfinite(a.laplacian) || !std::isfinite(a.tau) ||
                            std::abs(a.value-b.value) > 1e-10 ||
                            (a.gradient-b.gradient).cwiseAbs().maxCoeff() > 1e-10 ||
                            (a.hessian-b.hessian).cwiseAbs().maxCoeff() > 1e-10 ||
                            std::abs(a.laplacian-b.laplacian) > 1e-10 ||
                            std::abs(a.tau-b.tau) > 1e-10)
                            throw std::runtime_error("parallel derivative numerical mismatch");
                    }
                }
            }
        }
        std::cout << argv[file] << " parallel/serial derivatives OK\n";
    }
}
