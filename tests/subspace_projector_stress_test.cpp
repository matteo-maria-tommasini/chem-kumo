#include "chem/kumo/subspace_projector.hpp"
#include "chem/kumo/density.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

static double max_abs(const Eigen::MatrixXd& x) {
    return x.cwiseAbs().maxCoeff();
}

static void require(bool condition, const std::string& reason) {
    if (!condition) throw std::runtime_error(reason);
}

static void check(const Eigen::MatrixXd& s, Eigen::Index k, const std::string& name) {
    const Eigen::Index n = s.rows();
    if (k > n) return;
    // The first k coordinate vectors span a deterministic AO subspace.
    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(n, k);
    a.topRows(k).setIdentity();
    const Eigen::MatrixXd gram = a.transpose() * s * a;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(gram);
    require(eig.info() == Eigen::Success && eig.eigenvalues().minCoeff() > 0,
            "invalid Gram matrix");
    const double condition = eig.eigenvalues().maxCoeff() / eig.eigenvalues().minCoeff();

    SubspaceProjector projector(s, a);
    const Eigen::MatrixXd& q = projector.orthonormal_basis();
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(k, k);
    const double orth = max_abs(q.transpose() * s * q - identity);
    const double span = max_abs(projector.apply(a) - a);

    Eigen::MatrixXd x(n, 3);
    for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index j = 0; j < 3; ++j)
            x(i,j) = std::sin(0.013 * double(i+1) * double(j+2));
    const Eigen::MatrixXd projected = projector.apply(x);
    const double idem = max_abs(projector.apply(projected) - projected);
    const double residual_orth = max_abs(q.transpose() * s * (x - projected));

    // Nonorthogonal but well-conditioned change of basis within the subspace.
    Eigen::MatrixXd transform = identity;
    for (Eigen::Index i = 0; i + 1 < k; ++i) transform(i,i+1) = 0.05;
    SubspaceProjector changed(s, (a * transform).eval());
    const double invariant = max_abs(changed.apply(x) - projected);

    const double tolerance = 2e-7;
    require(orth < tolerance && span < tolerance && idem < tolerance &&
            residual_orth < tolerance && invariant < tolerance,
            "projector residual exceeds tolerance for " + name + " k=" + std::to_string(k));

    // Exactly dependent input must be rejected without silently dropping columns.
    if (k < n) {
        Eigen::MatrixXd dependent(n, k+1);
        dependent.leftCols(k) = a;
        dependent.col(k) = a.col(0);
        bool rejected = false;
        try { SubspaceProjector bad(s, dependent); (void)bad; }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "exact linear dependence accepted");

        // Nearly dependent input: make the last column a tiny independent perturbation.
        dependent.col(k) = a.col(0);
        dependent(k,k) = 1e-9;
        rejected = false;
        try { SubspaceProjector bad(s, dependent); (void)bad; }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "near-linear dependence accepted");
    }

    std::cout << std::scientific << std::setprecision(8)
              << name << " nAO=" << n << " k=" << k
              << " gram_condition=" << condition
              << " orth=" << orth << " span=" << span
              << " idem=" << idem << " residual_orth=" << residual_orth
              << " invariant=" << invariant << '\n';
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            FormattedCheckpoint fchk(argv[arg]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            const Eigen::MatrixXd s = evaluator.overlap_matrix();
            for (Eigen::Index k : {Eigen::Index{2}, Eigen::Index{10},
                                   Eigen::Index{50}, Eigen::Index{100}})
                check(s, k, argv[arg]);
        }
    } catch (const std::exception& e) {
        std::cerr << "projector stress: " << e.what() << '\n';
        return 1;
    }
}
