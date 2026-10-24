#pragma once
#include "chem/kumo/local_atomic_w.hpp"
#include "chem/kumo/density.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace chem::kumo {

struct NeighborWBlock {
    std::size_t center_a{}, center_b{};
    std::vector<Eigen::Index> ao_a, ao_b;
    std::vector<Eigen::MatrixXd> matrices;
    std::size_t evaluated_quartets{};
};

// Explicit first-neighbor graph; radius is in bohr, like FCHK coordinates.
[[nodiscard]] inline std::vector<std::pair<std::size_t,std::size_t>>
local_neighbor_pairs(const Molecule& molecule, double cutoff_bohr) {
    if (!std::isfinite(cutoff_bohr) || cutoff_bohr<=0)
        throw std::invalid_argument("neighbor W: invalid cutoff");
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    const auto& atoms=molecule.atoms();
    for (std::size_t a=0;a<atoms.size();++a)
        for (std::size_t b=a+1;b<atoms.size();++b)
            if ((atoms[a].position-atoms[b].position).norm()<cutoff_bohr)
                edges.emplace_back(a,b);
    return edges;
}

// Exact, unscreened AB block. The source density is molecular/global.
[[nodiscard]] inline NeighborWBlock local_neighbor_w_block(
    const BasisSet& basis, const std::vector<Eigen::MatrixXd>& densities,
    std::size_t center_a, std::size_t center_b) {
    const auto shells=local_shell_components(basis);
    std::vector<std::size_t> centers;
    for (const auto& shell:basis.shells()) {
        centers.push_back(shell.center);
        if (shell.is_sp) centers.push_back(shell.center);
    }
    if (centers.size()!=shells.size()) throw std::logic_error("neighbor W: shell mismatch");
    std::vector<Eigen::Index> offsets(shells.size()+1,0);
    for(std::size_t s=0;s<shells.size();++s)
        offsets[s+1]=offsets[s]+static_cast<Eigen::Index>(shells[s].size());
    const auto nao=offsets.back();
    for(const auto& p:densities) {
        if(p.rows()!=nao || p.cols()!=nao || !p.allFinite())
            throw std::invalid_argument("neighbor W: invalid density");
        if(nao && (p-p.transpose()).cwiseAbs().maxCoeff()>
           1e-10*std::max(1.0,p.cwiseAbs().maxCoeff()))
            throw std::invalid_argument("neighbor W: asymmetric density");
    }
    NeighborWBlock out;
    out.center_a=center_a;out.center_b=center_b;
    std::vector<std::size_t> sa,sb;
    for(std::size_t s=0;s<shells.size();++s) {
        if(centers[s]==center_a) {
            sa.push_back(s);
            for(auto i=offsets[s];i<offsets[s+1];++i) out.ao_a.push_back(i);
        }
        if(centers[s]==center_b) {
            sb.push_back(s);
            for(auto i=offsets[s];i<offsets[s+1];++i) out.ao_b.push_back(i);
        }
    }
    if(sa.empty()||sb.empty()) throw std::invalid_argument("neighbor W: unknown center");
    if(center_a==center_b) {
        const auto aa=local_atomic_w(basis,densities,center_a);
        out.matrices=aa.matrices;
        out.evaluated_quartets=aa.evaluated_quartets;
        return out;
    }
    const auto na=static_cast<Eigen::Index>(out.ao_a.size());
    const auto nb=static_cast<Eigen::Index>(out.ao_b.size());
    for(std::size_t ch=0;ch<densities.size();++ch)
        out.matrices.emplace_back(Eigen::MatrixXd::Zero(na,nb));
    std::vector<Eigen::Index> map_a(static_cast<std::size_t>(nao),-1),map_b(static_cast<std::size_t>(nao),-1);
    for(Eigen::Index i=0;i<na;++i)map_a[static_cast<std::size_t>(out.ao_a[static_cast<std::size_t>(i)])]=i;
    for(Eigen::Index j=0;j<nb;++j)map_b[static_cast<std::size_t>(out.ao_b[static_cast<std::size_t>(j)])]=j;
    for(const auto a:sa)for(const auto b:sb) {
        const auto target=prepare_local_shell_pair(shells[a],shells[b]);
        for(std::size_t c=0;c<shells.size();++c)for(std::size_t d=c;d<shells.size();++d) {
            const auto source=prepare_local_shell_pair(shells[c],shells[d]);
            const auto q=local_shell_quartet_paired(target,source);
            ++out.evaluated_quartets;
            for(std::size_t i=0;i<q.dimensions[0];++i)
            for(std::size_t j=0;j<q.dimensions[1];++j) {
                const auto li=map_a[static_cast<std::size_t>(offsets[a]+static_cast<Eigen::Index>(i))];
                const auto lj=map_b[static_cast<std::size_t>(offsets[b]+static_cast<Eigen::Index>(j))];
                for(std::size_t k=0;k<q.dimensions[2];++k)
                for(std::size_t l=0;l<q.dimensions[3];++l) {
                    const auto lambda=offsets[c]+static_cast<Eigen::Index>(k);
                    const auto sigma=offsets[d]+static_cast<Eigen::Index>(l);
                    if(lambda>sigma)continue;
                    const double v=(lambda==sigma?1.0:2.0)*q.at(i,j,k,l);
                    for(std::size_t ch=0;ch<densities.size();++ch)
                        out.matrices[ch](li,lj)+=v*densities[ch](lambda,sigma);
                }
            }
        }
    }
    return out;
}
} // namespace chem::kumo
