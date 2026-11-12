#include "chem/kumo/subspace_projector.hpp"
#include "chem/kumo/density.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

static void require(bool condition, const std::string& reason) {
    kumo_test::check(condition, reason);
}

static void check(const Eigen::MatrixXd& s, Eigen::Index k, const std::string& name) {
    const Eigen::Index n = s.rows();
    if (k > n) return;
    using kumo_test::max_abs_diff;
    const std::string where = " (" + name + " k=" + std::to_string(k) + ")";
    // The first k coordinate vectors span a deterministic AO subspace.
    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(n, k);
    a.topRows(k).setIdentity();
    const Eigen::MatrixXd gram = a.transpose() * s * a;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(gram);
    if (!kumo_test::check(eig.info() == Eigen::Success && eig.eigenvalues().minCoeff() > 0,
                          "Gram matrix of the subspace is positive definite" + where))
        return;
    const double condition = eig.eigenvalues().maxCoeff() / eig.eigenvalues().minCoeff();

    SubspaceProjector projector(s, a);
    const Eigen::MatrixXd& q = projector.orthonormal_basis();
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(k, k);
    const double orth = max_abs_diff(q.transpose() * s * q, identity);
    const double span = max_abs_diff(projector.apply(a), a);

    Eigen::MatrixXd x(n, 3);
    for (Eigen::Index i = 0; i < n; ++i)
        for (Eigen::Index j = 0; j < 3; ++j)
            x(i, j) = std::sin(0.013 * double(i + 1) * double(j + 2));
    const Eigen::MatrixXd projected = projector.apply(x);
    const double idem = max_abs_diff(projector.apply(projected), projected);
    const double residual_orth = (q.transpose() * s * (x - projected)).cwiseAbs().maxCoeff();

    // Nonorthogonal but well-conditioned change of basis within the subspace.
    Eigen::MatrixXd transform = identity;
    for (Eigen::Index i = 0; i + 1 < k; ++i)
        transform(i, i + 1) = 0.05;
    SubspaceProjector changed(s, (a * transform).eval());
    const double invariant = max_abs_diff(changed.apply(x), projected);

    const double tolerance = 2e-7;
    const auto below = [&](double error, const std::string& what) {
        require(error < tolerance, what + " within " + kumo_test::num(tolerance) +
                                       " (max_abs=" + kumo_test::num(error) + ")" + where);
    };
    below(orth, "basis is S-orthonormal");
    below(span, "projector preserves the subspace");
    below(idem, "projector is idempotent");
    below(residual_orth, "residual is S-orthogonal to the subspace");
    below(invariant, "projector is invariant under a change of basis");

    // Exactly dependent input must be rejected without silently dropping columns.
    if (k < n) {
        Eigen::MatrixXd dependent(n, k + 1);
        dependent.leftCols(k) = a;
        dependent.col(k) = a.col(0);
        bool rejected = false;
        try {
            SubspaceProjector bad(s, dependent);
            (void)bad;
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "exact linear dependence is rejected" + where);

        // Nearly dependent input: make the last column a tiny independent perturbation.
        dependent.col(k) = a.col(0);
        dependent(k, k) = 1e-9;
        rejected = false;
        try {
            SubspaceProjector bad(s, dependent);
            (void)bad;
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "near-linear dependence is rejected" + where);
    }

    std::cout << std::scientific << std::setprecision(8) << name << " nAO=" << n << " k=" << k
              << " gram_condition=" << condition << " orth=" << orth << " span=" << span
              << " idem=" << idem << " residual_orth=" << residual_orth
              << " invariant=" << invariant << '\n';
}

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        for (int arg = 1; arg < argc; ++arg) {
            FormattedCheckpoint fchk(argv[arg]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            const Eigen::MatrixXd s = evaluator.overlap_matrix();
            for (Eigen::Index k :
                 {Eigen::Index{2}, Eigen::Index{10}, Eigen::Index{50}, Eigen::Index{100}})
                check(s, k, argv[arg]);
        }
    });
}
