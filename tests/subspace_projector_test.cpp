#include "chem/kumo/subspace_projector.hpp"
#include "chem/kumo/density.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
static double max_abs(const Eigen::MatrixXd& x) { return x.cwiseAbs().maxCoeff(); }

static void validate(const Eigen::MatrixXd& s, const Eigen::MatrixXd& a,
                     double tolerance) {
    SubspaceProjector p(s, a);
    const Eigen::MatrixXd q = p.orthonormal_basis();
    const Eigen::MatrixXd pi = p.matrix();
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(a.cols(), a.cols());
    require(max_abs(q.transpose() * s * q - identity) < tolerance,
            "S-orthonormality failed");
    require(max_abs(pi * pi - pi) < tolerance, "idempotency failed");
    require(max_abs(pi.transpose() * s - s * pi) < tolerance,
            "metric self-adjointness failed");
    require(max_abs(p.apply(a) - a) < tolerance, "subspace preservation failed");

    Eigen::MatrixXd transform = Eigen::MatrixXd::Identity(a.cols(), a.cols());
    for (Eigen::Index i = 0; i < a.cols(); ++i)
        transform(i, i) = 0.8 + 0.1 * static_cast<double>(i);
    for (Eigen::Index i = 0; i + 1 < a.cols(); ++i)
        transform(i, i + 1) = 0.2;
    SubspaceProjector rotated(s, (a * transform).eval());
    require(max_abs(rotated.matrix() - pi) < tolerance,
            "change-of-basis invariance failed");

    Eigen::MatrixXd dependent(a.rows(), a.cols() + 1);
    dependent.leftCols(a.cols()) = a;
    dependent.col(a.cols()) = a.col(0);
    if (dependent.cols() <= dependent.rows()) {
        bool rejected = false;
        try { SubspaceProjector bad(s, dependent); (void)bad; }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "rank-deficient input accepted");
    }
}

int main(int argc, char** argv) {
    try {
        Eigen::MatrixXd s(4,4);
        s << 1.0,0.1,0.0,0.0,
             0.1,1.0,0.2,0.0,
             0.0,0.2,1.0,0.1,
             0.0,0.0,0.1,1.0;
        Eigen::MatrixXd a(4,2);
        a << 1.0,0.0, 0.2,0.3, 0.0,1.0, 0.1,0.2;
        validate(s,a,1e-11);
        if (argc < 2) return 2;
        for (int arg = 1; arg < argc; ++arg) {
            FormattedCheckpoint fchk(argv[arg]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            const auto metric = evaluator.overlap_matrix();
            // A small, deterministic subspace; do not assume any particular MO access API.
            Eigen::MatrixXd columns = Eigen::MatrixXd::Zero(metric.rows(), 2);
            columns(0,0) = 1.0;
            columns(1,1) = 1.0;
            validate(metric,columns,1e-9);
            std::cout << argv[arg] << " projector_dim=2 validated\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "subspace projector: " << e.what() << '\n';
        return 1;
    }
}
