#pragma once
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

struct LocalJointDiagonalization {
    Eigen::MatrixXd rotation;
    std::vector<Eigen::MatrixXd> transformed;
    double initial_offdiag = 0.0;
    double final_offdiag = 0.0;
    std::size_t sweeps = 0;
};

[[nodiscard]] inline double local_joint_offdiag(const std::vector<Eigen::MatrixXd>& matrices) {
    double result = 0.0;
    for (const Eigen::MatrixXd& g : matrices) {
        result += g.squaredNorm() - g.diagonal().squaredNorm();
    }
    return std::max(0.0, result);
}

// Real symmetric approximate joint diagonalization using orthogonal Jacobi rotations.
// The minimization is over the sum of squared off-diagonal entries of all matrices.
[[nodiscard]] inline LocalJointDiagonalization local_joint_diagonalize(
    const std::vector<Eigen::MatrixXd>& input,
    std::size_t max_sweeps = 200, double relative_tolerance = 1e-12) {
    if (input.empty()) throw std::invalid_argument("AJD: empty input");
    const Eigen::Index n = input.front().rows();
    if (n < 1 || !std::isfinite(relative_tolerance) || relative_tolerance <= 0)
        throw std::invalid_argument("AJD: invalid parameters");
    LocalJointDiagonalization out;
    out.rotation = Eigen::MatrixXd::Identity(n,n);
    for (const Eigen::MatrixXd& g : input) {
        if (g.rows() != n || g.cols() != n || !g.allFinite() ||
            (g-g.transpose()).norm() > 1e-10 * std::max(1.0,g.norm()))
            throw std::invalid_argument("AJD: matrices must be finite and symmetric");
        out.transformed.push_back((0.5*(g+g.transpose())).eval());
    }
    out.initial_offdiag = local_joint_offdiag(out.transformed);
    double previous = out.initial_offdiag;
    constexpr double pi = 3.14159265358979323846;
    for (std::size_t sweep=0; sweep<max_sweeps; ++sweep) {
        bool changed = false;
        for (Eigen::Index p=0;p<n;++p) for(Eigen::Index q=p+1;q<n;++q) {
            double aa=0.0, bb=0.0, ab=0.0;
            for (const Eigen::MatrixXd& g : out.transformed) {
                const double a=g(p,q);
                const double b=0.5*(g(p,p)-g(q,q));
                aa+=a*a; bb+=b*b; ab+=a*b;
            }
            // Off-diagonal after rotation: a*cos(2 theta)+b*sin(2 theta).
            const double angle=0.25*std::atan2(2.0*ab,aa-bb);
            const auto objective=[&](double theta) {
                const double c=std::cos(2*theta),s=std::sin(2*theta);
                return aa*c*c+2*ab*c*s+bb*s*s;
            };
            const double old=aa;
            const double candidate1=angle;
            const double candidate2=angle+pi/4.0;
            const double theta=objective(candidate1)<objective(candidate2)?candidate1:candidate2;
            const double best=objective(theta);
            if (old-best <= 1e-15*std::max({1.0,aa,bb}) ||
                old-best <= relative_tolerance*std::max(1.0,out.initial_offdiag)) continue;
            const double c=std::cos(theta),s=std::sin(theta);
            Eigen::MatrixXd r=Eigen::MatrixXd::Identity(n,n);
            r(p,p)=c; r(q,q)=c; r(p,q)=s; r(q,p)=-s;
            for (Eigen::MatrixXd& g:out.transformed) {
                const Eigen::MatrixXd updated=(r.transpose()*g*r).eval();
                g=0.5*(updated+updated.transpose());
            }
            out.rotation=(out.rotation*r).eval();
            changed=true;
        }
        out.sweeps=sweep+1;
        const double current=local_joint_offdiag(out.transformed);
        if (!changed || previous-current <= relative_tolerance*std::max(1.0,out.initial_offdiag)) break;
        previous=current;
    }
    out.final_offdiag=local_joint_offdiag(out.transformed);
    return out;
}
}
