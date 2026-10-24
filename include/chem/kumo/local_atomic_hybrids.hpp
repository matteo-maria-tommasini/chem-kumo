#pragma once

#include "chem/kumo/local_atomic_w.hpp"
#include "chem/kumo/basis_evaluator.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chem::kumo {

struct LocalAtomicHybridChannel {
    Eigen::VectorXd eigenvalues; // descending weighted-density eigenvalues
    Eigen::MatrixXd coefficients; // nlocal x rank, C^T S C = I
    double residual_max = 0.0; // ||W C - S C Lambda||_max
    double orthogonality_max = 0.0;
};

struct LocalAtomicHybrids {
    std::size_t center = 0;
    std::vector<Eigen::Index> ao_indices;
    Eigen::MatrixXd overlap;
    Eigen::Index rank = 0;
    Eigen::VectorXd overlap_eigenvalues;
    std::vector<LocalAtomicHybridChannel> channels;
};

// The supplied overlap is the GLOBAL analytic AO overlap in the same AO order
// as local_atomic_w(). Only the on-center principal submatrix is used.
// Eigenvectors within exactly/near-degenerate eigenspaces are not unique.
[[nodiscard]] inline LocalAtomicHybrids local_atomic_hybrids(
    const LocalAtomicW& atomic, const Eigen::MatrixXd& global_overlap,
    double relative_rank_threshold = 1e-10) {
    if (!(relative_rank_threshold > 0.0) || !(relative_rank_threshold < 1.0) ||
        !std::isfinite(relative_rank_threshold))
        throw std::invalid_argument("atomic hybrids: invalid rank threshold");
    if (global_overlap.rows() != global_overlap.cols() || !global_overlap.allFinite())
        throw std::invalid_argument("atomic hybrids: invalid global overlap");
    const Eigen::Index n = static_cast<Eigen::Index>(atomic.ao_indices.size());
    if (n == 0) throw std::invalid_argument("atomic hybrids: empty AO selection");
    LocalAtomicHybrids result;
    result.center = atomic.center;
    result.ao_indices = atomic.ao_indices;
    result.overlap.resize(n,n);
    std::vector<bool> seen(static_cast<std::size_t>(global_overlap.rows()),false);
    for (const auto index : atomic.ao_indices)
        if (index<0 || index>=global_overlap.rows())
            throw std::invalid_argument("atomic hybrids: AO index out of range");
    for (Eigen::Index i=0;i<n;++i) {
        const auto gi=atomic.ao_indices[static_cast<std::size_t>(i)];
        if (gi<0 || gi>=global_overlap.rows() || seen[static_cast<std::size_t>(gi)])
            throw std::invalid_argument("atomic hybrids: invalid or duplicate AO index");
        seen[static_cast<std::size_t>(gi)]=true;
        for (Eigen::Index j=0;j<n;++j)
            result.overlap(i,j)=global_overlap(gi,atomic.ao_indices[static_cast<std::size_t>(j)]);
    }
    const auto& s=result.overlap;
    const double sscale=std::max(1.0,s.cwiseAbs().maxCoeff());
    if ((s-s.transpose()).cwiseAbs().maxCoeff()>1e-10*sscale)
        throw std::invalid_argument("atomic hybrids: asymmetric overlap");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> metric(0.5*(s+s.transpose()));
    if(metric.info()!=Eigen::Success) throw std::runtime_error("atomic hybrids: overlap eigensolver failed");
    result.overlap_eigenvalues=metric.eigenvalues();
    const double largest=metric.eigenvalues()[n-1];
    if(!(largest>0.0)) throw std::invalid_argument("atomic hybrids: nonpositive overlap");
    if(metric.eigenvalues()[0] < -1e-10*largest)
        throw std::invalid_argument("atomic hybrids: indefinite overlap");
    const double cutoff=relative_rank_threshold*largest;
    for(Eigen::Index i=0;i<n;++i) if(metric.eigenvalues()[i]>cutoff) ++result.rank;
    if(result.rank==0) throw std::invalid_argument("atomic hybrids: zero numerical rank");
    Eigen::MatrixXd x(n,result.rank);
    Eigen::Index col=0;
    for(Eigen::Index i=0;i<n;++i) if(metric.eigenvalues()[i]>cutoff)
        x.col(col++)=metric.eigenvectors().col(i)/std::sqrt(metric.eigenvalues()[i]);
    for(const auto& raw:atomic.matrices) {
        if(raw.rows()!=n || raw.cols()!=n || !raw.allFinite())
            throw std::invalid_argument("atomic hybrids: invalid W block");
        const double wscale=std::max(1.0,raw.cwiseAbs().maxCoeff());
        if((raw-raw.transpose()).cwiseAbs().maxCoeff()>1e-10*wscale)
            throw std::invalid_argument("atomic hybrids: asymmetric W block");
        const Eigen::MatrixXd w=0.5*(raw+raw.transpose());
        const Eigen::MatrixXd transformed=x.transpose()*w*x;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(transformed);
        if(eig.info()!=Eigen::Success) throw std::runtime_error("atomic hybrids: W eigensolver failed");
        LocalAtomicHybridChannel ch;
        ch.eigenvalues.resize(result.rank);
        ch.coefficients.resize(n,result.rank);
        for(Eigen::Index i=0;i<result.rank;++i) {
            const Eigen::Index j=result.rank-1-i;
            ch.eigenvalues[i]=eig.eigenvalues()[j];
            ch.coefficients.col(i)=x*eig.eigenvectors().col(j);
            Eigen::Index pivot=0;
            ch.coefficients.col(i).cwiseAbs().maxCoeff(&pivot);
            if(ch.coefficients(pivot,i)<0.0) ch.coefficients.col(i)*=-1.0;
        }
        const Eigen::MatrixXd gram=ch.coefficients.transpose()*s*ch.coefficients;
        ch.orthogonality_max=(gram-Eigen::MatrixXd::Identity(result.rank,result.rank)).cwiseAbs().maxCoeff();
        // For rank-truncated S, the full-space residual need not vanish:
        // report the Galerkin residual in the retained orthonormal subspace.
        const Eigen::MatrixXd residual=x.transpose()*(w*ch.coefficients-s*ch.coefficients*ch.eigenvalues.asDiagonal());
        ch.residual_max=residual.cwiseAbs().maxCoeff();
        result.channels.push_back(std::move(ch));
    }
    return result;
}

[[nodiscard]] inline LocalAtomicHybrids local_atomic_hybrids(
    const BasisSet& basis, const LocalAtomicW& atomic,
    double relative_rank_threshold=1e-10) {
    return local_atomic_hybrids(atomic,BasisEvaluator(basis).overlap_matrix(),relative_rank_threshold);
}

} // namespace chem::kumo
