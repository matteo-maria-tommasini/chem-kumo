#include "chem/kumo/basis_evaluator.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int a = 1; a < argc; ++a) {
            chem::kumo::FormattedCheckpoint fchk(argv[a]);
            chem::kumo::BasisEvaluator evaluator{chem::kumo::BasisSet(fchk)};
            const Eigen::MatrixXd S = evaluator.overlap_matrix();
            if (!S.allFinite() || (S-S.transpose()).cwiseAbs().maxCoeff() > 1e-12)
                throw std::runtime_error("overlap not finite or symmetric");
            const double symmetry = (S - S.transpose()).cwiseAbs().maxCoeff();
            const double diagonal_error = (S.diagonal().array() - 1.0).abs().maxCoeff();
            // Full symmetric eigendecomposition is intentionally limited to small bases.
            // This diagnostic must remain usable for 1000+ AO regression systems.
            if (S.rows() <= 128) {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S, Eigen::EigenvaluesOnly);
                if (solver.info() != Eigen::Success)
                    throw std::runtime_error("overlap eigensolver failed");
                const double lambda_min = solver.eigenvalues()[0];
                const double lambda_max = solver.eigenvalues()[solver.eigenvalues().size()-1];
                std::cout << std::scientific << std::setprecision(12)
                          << argv[a] << " overlap_symmetry=" << symmetry
                          << " ao_diagonal_error=" << diagonal_error
                          << " lambda_min=" << lambda_min
                          << " lambda_max=" << lambda_max
                          << " condition=" << (lambda_min > 0.0 ? lambda_max/lambda_min : std::numeric_limits<double>::infinity())
                          << '\n';
                if (!(lambda_min > 0.0)) throw std::runtime_error("overlap not positive definite");
            } else {
                std::cout << std::scientific << std::setprecision(12)
                          << argv[a] << " overlap_symmetry=" << symmetry
                          << " ao_diagonal_error=" << diagonal_error
                          << " eigenvalues=skipped_large_basis" << '\n';
            }
            const auto test_spin = [&](chem::kumo::OrbitalSpin spin) {
                chem::kumo::MolecularOrbitals orbitals(fchk, spin);
                const auto& C = orbitals.coefficients();
                const Eigen::MatrixXd metric = C.transpose() * S * C;
                const double error = (metric - Eigen::MatrixXd::Identity(metric.rows(), metric.cols())).cwiseAbs().maxCoeff();
                const Eigen::VectorXd diagonal = metric.diagonal();
                const double diag_error = (diagonal.array() - 1.0).abs().maxCoeff();
                Eigen::MatrixXd offdiag = metric;
                offdiag.diagonal().setZero();
                const double offdiag_error = offdiag.cwiseAbs().maxCoeff();
                std::cout << std::scientific << std::setprecision(12)
                          << argv[a] << " spin=" << (spin == chem::kumo::OrbitalSpin::Alpha ? "alpha" : "beta")
                          << " nAO=" << S.rows() << " orthonormality_max=" << error
                          << " diagonal_error=" << diag_error
                          << " offdiagonal_error=" << offdiag_error << '\n';
                if (!std::isfinite(error) || error > 1e-5)
                    throw std::runtime_error("MO orthonormality regression");
            };
            test_spin(chem::kumo::OrbitalSpin::Alpha);
            // Restricted FCHK may omit beta coefficients.
            if (fchk.contains("Beta MO coefficients")) test_spin(chem::kumo::OrbitalSpin::Beta);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
