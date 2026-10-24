#pragma once

#include "chem/kumo/local_density_overlap_multichannel.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

struct LocalAtomicW {
    std::size_t center = 0;
    std::vector<Eigen::Index> ao_indices;
    std::vector<Eigen::MatrixXd> matrices;
    std::size_t evaluated_quartets = 0;
    double total_seconds = 0.0;
};

// Unscreened reference for the on-center AO block, using the GLOBAL density.
// Only target shell pairs on the requested center are evaluated. The source
// shell pairs cover the complete molecular basis. No global W is allocated.
[[nodiscard]] inline LocalAtomicW local_atomic_w(
    const BasisSet& basis, const std::vector<Eigen::MatrixXd>& densities,
    std::size_t center) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    const auto shells = local_shell_components(basis);
    std::vector<std::size_t> centers;
    centers.reserve(shells.size());
    for (const auto& shell : basis.shells()) {
        centers.push_back(shell.center);
        if (shell.is_sp) centers.push_back(shell.center);
    }
    if (centers.size() != shells.size())
        throw std::logic_error("atomic W: shell/component center mismatch");
    std::vector<Eigen::Index> offsets(shells.size()+1, 0);
    for (std::size_t i=0; i<shells.size(); ++i)
        offsets[i+1] = offsets[i] + static_cast<Eigen::Index>(shells[i].size());
    const Eigen::Index nao = offsets.back();
    if (center >= basis.shells().size() &&
        std::find(centers.begin(),centers.end(),center)==centers.end())
        throw std::invalid_argument("atomic W: unknown center");
    for (const auto& p : densities) {
        if (p.rows()!=nao || p.cols()!=nao || !p.allFinite())
            throw std::invalid_argument("atomic W: invalid density");
        if (nao && (p-p.transpose()).cwiseAbs().maxCoeff() >
                1e-10*std::max(1.0,p.cwiseAbs().maxCoeff()))
            throw std::invalid_argument("atomic W: asymmetric density");
    }
    LocalAtomicW out;
    out.center=center;
    std::vector<std::size_t> target_shells;
    for (std::size_t s=0;s<shells.size();++s) if (centers[s]==center) {
        target_shells.push_back(s);
        for (Eigen::Index i=offsets[s];i<offsets[s+1];++i)
            out.ao_indices.push_back(i);
    }
    if (target_shells.empty()) throw std::invalid_argument("atomic W: unknown center");
    const Eigen::Index nloc=static_cast<Eigen::Index>(out.ao_indices.size());
    for (std::size_t k=0;k<densities.size();++k)
        out.matrices.emplace_back(Eigen::MatrixXd::Zero(nloc,nloc));
    std::vector<Eigen::Index> local_of_global(static_cast<std::size_t>(nao),-1);
    for (Eigen::Index i=0;i<nloc;++i)
        local_of_global[static_cast<std::size_t>(out.ao_indices[static_cast<std::size_t>(i)])]=i;
    struct Pair { std::size_t a,b; };
    std::vector<Pair> source_pairs;
    for (std::size_t c=0;c<shells.size();++c)
        for (std::size_t d=c;d<shells.size();++d)
            source_pairs.push_back({c,d});
    for (std::size_t ti=0;ti<target_shells.size();++ti)
        for (std::size_t tj=ti;tj<target_shells.size();++tj) {
            const auto a=target_shells[ti], b=target_shells[tj];
            const auto target=prepare_local_shell_pair(shells[a],shells[b]);
            for (const auto [c,d] : source_pairs) {
                const auto source=prepare_local_shell_pair(shells[c],shells[d]);
                const auto q=local_shell_quartet_paired(target,source);
                ++out.evaluated_quartets;
                                for (std::size_t i=0;i<q.dimensions[0];++i)
                    for (std::size_t j=0;j<q.dimensions[1];++j) {
                        const auto mu=offsets[a]+static_cast<Eigen::Index>(i);
                        const auto nu=offsets[b]+static_cast<Eigen::Index>(j);
                        if (mu>nu) continue;
                        const auto li=local_of_global[static_cast<std::size_t>(mu)];
                        const auto lj=local_of_global[static_cast<std::size_t>(nu)];
                        for (std::size_t k=0;k<q.dimensions[2];++k)
                            for (std::size_t l=0;l<q.dimensions[3];++l) {
                                const auto lambda=offsets[c]+static_cast<Eigen::Index>(k);
                                const auto sigma=offsets[d]+static_cast<Eigen::Index>(l);
                                if (lambda>sigma) continue;
                                const double v=((lambda==sigma)?1.0:2.0)*q.at(i,j,k,l);
                                for (std::size_t ch=0;ch<densities.size();++ch)
                                    out.matrices[ch](li,lj)+=v*densities[ch](lambda,sigma);
                            }
                    }
            }
        }
    for (auto& w : out.matrices)
        w=Eigen::MatrixXd(w.selfadjointView<Eigen::Upper>());
    out.total_seconds=std::chrono::duration<double>(Clock::now()-start).count();
    return out;
}

} // namespace chem::kumo
