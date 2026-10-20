#include "chem/kumo/weighted_overlap.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

static double error(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
    return (a-b).cwiseAbs().maxCoeff();
}
static void require(bool ok, const char* why) {
    if (!ok) throw std::runtime_error(why);
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            FormattedCheckpoint fchk(argv[arg]);
            BasisEvaluator evaluator{BasisSet(fchk)};
            constexpr Eigen::Index n = 17;
            Eigen::Matrix<double,3,Eigen::Dynamic> points(3,n);
            Eigen::VectorXd weights(n), field(n), ones = Eigen::VectorXd::Ones(n);
            for (Eigen::Index j=0; j<n; ++j) {
                points.col(j) = Eigen::Vector3d(0.13*j-0.8, 0.17*j-0.5, 0.09*j+0.2);
                weights[j] = 0.01 + 0.002*j;
                field[j] = 0.3 + 0.07*j;
            }
            Eigen::MatrixXd ao(evaluator.basis_size(), n);
            evaluator.evaluate_batch(points, ao);
            const Eigen::VectorXd effective = (weights.array()*field.array()).matrix();
            const Eigen::MatrixXd reference = ao * effective.asDiagonal() * ao.transpose();
            const Eigen::MatrixXd blocked = weighted_overlap_matrix(evaluator,points,weights,field,5);
            const Eigen::MatrixXd singleton = weighted_overlap_matrix(evaluator,points,weights,field,1);
            const double tolerance = 3e-12 * std::max(1.0,reference.cwiseAbs().maxCoeff());
            require(error(reference,blocked) < tolerance, "blocked reference mismatch");
            require(error(reference,singleton) < tolerance, "single-point reference mismatch");
            require(error(blocked,blocked.transpose()) < tolerance, "symmetry mismatch");
            const Eigen::VectorXd double_field = 2.0 * field;
            const Eigen::MatrixXd twice = weighted_overlap_matrix(evaluator,points,weights,double_field,7);
            require(error(twice,2.0*reference) < 2*tolerance, "linearity mismatch");
            const Eigen::VectorXd zero_field = Eigen::VectorXd::Zero(n);
            const Eigen::MatrixXd zero = weighted_overlap_matrix(evaluator,points,weights,zero_field,7);
            require(zero.cwiseAbs().maxCoeff() == 0.0, "zero field mismatch");
            // Positive semidefiniteness without a costly full eigendecomposition.
            for (int k=0; k<4; ++k) {
                Eigen::VectorXd v = Eigen::VectorXd::Zero(evaluator.basis_size());
                for (Eigen::Index i=0; i<v.size(); ++i) v[i] = std::sin(0.17*(i+1)*(k+1));
                require(v.dot(blocked*v) >= -tolerance*v.squaredNorm(), "PSD mismatch");
            }
            // Signed weights are permitted and must not be silently clipped.
            const Eigen::VectorXd negative_field = -field;
            const Eigen::MatrixXd negative = weighted_overlap_matrix(evaluator,points,weights,negative_field,6);
            require(error(negative,-reference) < tolerance, "signed field mismatch");
            WeightedOverlapAccumulator accumulator(evaluator);
            accumulator.add(points,weights,field);
            require(error(accumulator.matrix(),reference) < tolerance, "accumulator mismatch");
            accumulator.reset();
            require(accumulator.matrix().cwiseAbs().maxCoeff() == 0.0, "reset mismatch");
            bool rejected = false;
            try { (void)weighted_overlap_matrix(evaluator,points,weights,field,0); }
            catch (const std::invalid_argument&) { rejected=true; }
            require(rejected,"zero block accepted");
            std::cout << argv[arg] << " nAO=" << evaluator.basis_size()
                      << " points=" << n << " max_abs=" << error(reference,blocked) << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << "weighted overlap: " << e.what() << '\n'; return 1;
    }
}
