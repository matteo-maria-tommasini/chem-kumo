#include "chem/kumo/subspace_projector.hpp"
#include "chem/kumo/density.hpp"
#include "test_support.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

static void require(bool condition, const std::string& message) {
    kumo_test::check(condition, message);
}

// Pass if max |actual - expected| < tolerance; the message shows both.
static void require_below(double error, double tolerance, const std::string& what) {
    kumo_test::check(error < tolerance, what + " (max_abs=" + kumo_test::num(error) +
                                            ", tolerance " + kumo_test::num(tolerance) + ")");
}

static void validate(const Eigen::MatrixXd& s, const Eigen::MatrixXd& a, double tolerance) {
    using kumo_test::max_abs_diff;
    SubspaceProjector p(s, a);
    const Eigen::MatrixXd q = p.orthonormal_basis();
    const Eigen::MatrixXd pi = p.matrix();
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(a.cols(), a.cols());
    require_below(max_abs_diff(q.transpose() * s * q, identity), tolerance,
                  "basis is S-orthonormal");
    require_below(max_abs_diff(pi * pi, pi), tolerance, "projector is idempotent");
    require_below(max_abs_diff(pi.transpose() * s, s * pi), tolerance,
                  "projector is self-adjoint in the metric");
    require_below(max_abs_diff(p.apply(a), a), tolerance, "projector preserves the subspace");

    Eigen::MatrixXd transform = Eigen::MatrixXd::Identity(a.cols(), a.cols());
    for (Eigen::Index i = 0; i < a.cols(); ++i)
        transform(i, i) = 0.8 + 0.1 * static_cast<double>(i);
    for (Eigen::Index i = 0; i + 1 < a.cols(); ++i)
        transform(i, i + 1) = 0.2;
    SubspaceProjector rotated(s, (a * transform).eval());
    require_below(max_abs_diff(rotated.matrix(), pi), tolerance,
                  "projector is invariant under a change of basis");

    Eigen::MatrixXd dependent(a.rows(), a.cols() + 1);
    dependent.leftCols(a.cols()) = a;
    dependent.col(a.cols()) = a.col(0);
    if (dependent.cols() <= dependent.rows()) {
        bool rejected = false;
        try {
            SubspaceProjector bad(s, dependent);
            (void)bad;
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "rank-deficient input is rejected");
    }
}

static void validate_invalid_inputs() {
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(4, 4);
    const Eigen::MatrixXd columns = identity.leftCols(2);

    auto rejected = [](const Eigen::MatrixXd& metric, const Eigen::MatrixXd& subspace) {
        try {
            const SubspaceProjector projector(metric, subspace);
            (void)projector;
        } catch (const std::invalid_argument&) {
            return true;
        }
        return false;
    };

    Eigen::MatrixXd nonsymmetric = identity;
    nonsymmetric(0, 1) = 0.1;
    require(rejected(nonsymmetric, columns), "a nonsymmetric metric is rejected");

    Eigen::MatrixXd indefinite = identity;
    indefinite(0, 0) = -1.0;
    require(rejected(indefinite, columns), "an indefinite metric is rejected");

    Eigen::MatrixXd singular = identity;
    singular(0, 0) = 0.0;
    require(rejected(singular, columns), "a singular metric is rejected");

    Eigen::MatrixXd nonfinite_metric = identity;
    nonfinite_metric(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require(rejected(nonfinite_metric, columns), "a nonfinite metric is rejected");

    Eigen::MatrixXd nonfinite_columns = columns;
    nonfinite_columns(0, 0) = std::numeric_limits<double>::infinity();
    require(rejected(identity, nonfinite_columns), "nonfinite columns are rejected");

    const SubspaceProjector projector(identity, columns);
    Eigen::MatrixXd vectors = Eigen::MatrixXd::Ones(4, 2);
    vectors(0, 0) = std::numeric_limits<double>::quiet_NaN();
    bool failed = false;
    try {
        (void)projector.apply(vectors);
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    require(failed, "a NaN projection input is rejected");

    vectors(0, 0) = std::numeric_limits<double>::infinity();
    failed = false;
    try {
        (void)projector.apply(vectors);
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    require(failed, "an infinite projection input is rejected");

    failed = false;
    try {
        (void)projector.apply(Eigen::MatrixXd::Ones(3, 2));
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    require(failed, "wrong projection dimensions are rejected");

    failed = false;
    try {
        (void)SubspaceProjector(identity, columns, 0.0);
    } catch (const std::invalid_argument&) {
        failed = true;
    }
    require(failed, "an invalid rank tolerance is rejected");
}

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        Eigen::MatrixXd s(4, 4);
        s << 1.0, 0.1, 0.0, 0.0, 0.1, 1.0, 0.2, 0.0, 0.0, 0.2, 1.0, 0.1, 0.0, 0.0, 0.1, 1.0;
        Eigen::MatrixXd a(4, 2);
        a << 1.0, 0.0, 0.2, 0.3, 0.0, 1.0, 0.1, 0.2;
        validate(s, a, 1e-11);
        validate_invalid_inputs();
        for (int arg = 1; arg < argc; ++arg) {
            FormattedCheckpoint fchk(argv[arg]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            const auto metric = evaluator.overlap_matrix();
            // A small, deterministic subspace; do not assume any particular MO access API.
            Eigen::MatrixXd columns = Eigen::MatrixXd::Zero(metric.rows(), 2);
            columns(0, 0) = 1.0;
            columns(1, 1) = 1.0;
            validate(metric, columns, 1e-9);
            std::cout << argv[arg] << " projector_dim=2 validated\n";
        }
    });
}
