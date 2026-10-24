#pragma once
#include "chem/kumo/basis.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/density.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace chem::kumo {
struct DensityPairScore {
    std::size_t a{}, b{};
    double distance_bohr{}, raw_frobenius{}, metric_frobenius{}, normalized{};
};

// Atom-local AO indices follow the original FCHK shell order, including SP shells.
[[nodiscard]] inline std::vector<std::vector<Eigen::Index>> local_atom_ao_indices(
    const BasisSet& basis, std::size_t natoms) {
    std::vector<std::vector<Eigen::Index>> groups(natoms);
    Eigen::Index offset=0;
    for (const auto& shell : basis.shells()) {
        if (shell.center>=natoms) throw std::invalid_argument("pair graph: shell center out of range");
        for (std::size_t i=0;i<shell.function_count();++i)
            groups[shell.center].push_back(offset++);
    }
    if (offset!=static_cast<Eigen::Index>(basis.function_count()))
        throw std::logic_error("pair graph: AO ordering mismatch");
    return groups;
}

[[nodiscard]] inline Eigen::MatrixXd local_ao_submatrix(
    const Eigen::MatrixXd& m, const std::vector<Eigen::Index>& a,
    const std::vector<Eigen::Index>& b) {
    Eigen::MatrixXd out(a.size(),b.size());
    for (Eigen::Index i=0;i<out.rows();++i)
        for (Eigen::Index j=0;j<out.cols();++j)
            out(i,j)=m(a[static_cast<std::size_t>(i)],b[static_cast<std::size_t>(j)]);
    return out;
}

// S^(1/2) uses a positive-definite local metric; no silent rank truncation.
[[nodiscard]] inline Eigen::MatrixXd local_metric_sqrt(const Eigen::MatrixXd& s) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
    if (es.info()!=Eigen::Success || !es.eigenvalues().allFinite())
        throw std::runtime_error("pair graph: overlap eigensolver failed");
    const double scale=std::max(1.0,es.eigenvalues().maxCoeff());
    if (es.eigenvalues().minCoeff()<=1e-12*scale)
        throw std::runtime_error("pair graph: local overlap near singular; explicit rank policy required");
    return es.eigenvectors()*es.eigenvalues().cwiseSqrt().asDiagonal()*es.eigenvectors().transpose();
}

[[nodiscard]] inline std::vector<DensityPairScore> local_density_pair_graph(
    const BasisSet& basis, const Molecule& molecule, const Eigen::MatrixXd& density,
    const Eigen::MatrixXd& overlap) {
    const auto nao=static_cast<Eigen::Index>(basis.function_count());
    if (density.rows()!=nao || density.cols()!=nao || overlap.rows()!=nao || overlap.cols()!=nao ||
        !density.allFinite() || !overlap.allFinite())
        throw std::invalid_argument("pair graph: invalid matrix dimensions or nonfinite values");
    const auto symcheck=[](const Eigen::MatrixXd& m) {
        return (m-m.transpose()).cwiseAbs().maxCoeff() <=
            1e-10*std::max(1.0,m.cwiseAbs().maxCoeff());
    };
    if (!symcheck(density)||!symcheck(overlap))
        throw std::invalid_argument("pair graph: asymmetric matrix");
    const auto groups=local_atom_ao_indices(basis,molecule.atoms().size());
    std::vector<Eigen::MatrixXd> roots;
    std::vector<double> diagonal_norms;
    for (const auto& g:groups) {
        if (g.empty()) throw std::invalid_argument("pair graph: atom without AO");
        auto r=local_metric_sqrt(local_ao_submatrix(overlap,g,g));
        const Eigen::MatrixXd p_aa =
            local_ao_submatrix(density, g, g);
        const Eigen::MatrixXd d =
            (r * p_aa * r).eval();
        diagonal_norms.push_back(d.norm());
        roots.push_back(std::move(r));
    }
    std::vector<DensityPairScore> result;
    for (std::size_t a=0;a<groups.size();++a)
        for (std::size_t b=a+1;b<groups.size();++b) {
            const auto p=local_ao_submatrix(density,groups[a],groups[b]);
//            const double eta=(roots[a]*p*roots[b]).norm();
const Eigen::MatrixXd intermediate =
    (roots[a] * p).eval();
const Eigen::MatrixXd metric_block =
    (intermediate * roots[b]).eval();
const double eta = metric_block.norm();
            const double denom=std::sqrt(diagonal_norms[a]*diagonal_norms[b]);
            const double gamma=denom>0 ? eta/denom : 0.0;
            result.push_back({a,b,(molecule.atoms()[a].position-molecule.atoms()[b].position).norm(),
                p.norm(),eta,gamma});
        }
    return result;
}
} // namespace chem::kumo
