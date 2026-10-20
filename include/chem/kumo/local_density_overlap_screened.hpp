#pragma once

#include "chem/kumo/local_density_quartet_contraction.hpp"

#include "chem/kumo/local_density_overlap_symmetric_paired.hpp"
#include "chem/kumo/local_pair_bounds.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chem::kumo {

struct LocalScreeningStatistics {
    std::size_t total_quartets = 0;
    std::size_t skipped_quartets = 0;
    double max_error_bound = 0.0;
    double norm_seconds = 0.0;
    double selection_seconds = 0.0;
    double pair_preparation_seconds = 0.0;
    double contraction_seconds = 0.0;
};

struct LocalScreenedDensityOverlap {
    Eigen::MatrixXd matrix;
    LocalScreeningStatistics statistics;
};

/// Density-aware CS screening with an absolute per-entry error budget.
/// The bound is rigorous in exact arithmetic; floating-point rounding is not certified.
/// No global four-index tensor is constructed.
[[nodiscard]] inline LocalScreenedDensityOverlap local_density_overlap_screened(
    const BasisSet& basis, const Eigen::Ref<const Eigen::MatrixXd>& density,
    double epsilon) {
    using Clock = std::chrono::steady_clock;
    if (!std::isfinite(epsilon) || epsilon < 0.0)
        throw std::invalid_argument("local_density_overlap_screened: invalid error budget");
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a+1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao = offsets.back();
    if (density.rows() != nao || density.cols() != nao || !density.allFinite())
        throw std::invalid_argument("local_density_overlap_screened: invalid density");
    LocalScreenedDensityOverlap out;
    out.matrix = Eigen::MatrixXd::Zero(nao, nao);
    if (!nao) return out;
    const double scale = density.cwiseAbs().maxCoeff();
    if ((density-density.transpose()).cwiseAbs().maxCoeff() > 1e-10*std::max(1.0,scale))
        throw std::invalid_argument("local_density_overlap_screened: asymmetric density");

    struct Pair { std::size_t a,b; double bmax=0, dsum=0; };
    std::vector<Pair> pairs;
    for (std::size_t a=0; a<shells.size(); ++a)
        for (std::size_t b=a; b<shells.size(); ++b)
            pairs.push_back({a,b});
    const std::size_t np=pairs.size();
    out.statistics.total_quartets=np*(np+1)/2;
    if (scale == 0.0) {
        out.statistics.skipped_quartets=out.statistics.total_quartets;
        return out;
    }

    const auto t0=Clock::now();
    const Eigen::MatrixXd norms=local_ao_product_norms(basis);
    const auto t1=Clock::now();
    out.statistics.norm_seconds=std::chrono::duration<double>(t1-t0).count();
    for (auto& pair:pairs) {
        for (Eigen::Index mu=offsets[pair.a];mu<offsets[pair.a+1];++mu)
            for (Eigen::Index nu=offsets[pair.b];nu<offsets[pair.b+1];++nu) {
                const double b=norms(mu,nu);
                pair.bmax=std::max(pair.bmax,b);
                pair.dsum+=std::abs(density(mu,nu))*b;
            }
    }
    struct Candidate { std::size_t i,j; long double to_i,to_j; };
    std::vector<Candidate> candidates;
    candidates.reserve(out.statistics.total_quartets);
    for (std::size_t i=0;i<np;++i)
        for (std::size_t j=i;j<np;++j) {
            const auto& x=pairs[i];
            const auto& y=pairs[j];
            const double mi=x.a==x.b?1.0:2.0;
            const double mj=y.a==y.b?1.0:2.0;
            candidates.push_back({i,j,
                static_cast<long double>(mj)*x.bmax*y.dsum,
                i==j?0.0L:static_cast<long double>(mi)*y.bmax*x.dsum});
        }
    std::vector<std::size_t> order(candidates.size());
    std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](std::size_t ai,std::size_t bi){
        const auto& a=candidates[ai];
        const auto& b=candidates[bi];
        const auto va=std::max(a.to_i,a.to_j);
        const auto vb=std::max(b.to_i,b.to_j);
        return va==vb?ai<bi:va<vb;
    });
    std::vector<long double> spent(np,0.0L);
    std::vector<unsigned char> skip(candidates.size(),0);
    for (const std::size_t index:order) {
        const auto& c=candidates[index];
        if (spent[c.i]+c.to_i>static_cast<long double>(epsilon)) continue;
        if (c.i!=c.j && spent[c.j]+c.to_j>static_cast<long double>(epsilon)) continue;
        spent[c.i]+=c.to_i;
        if (c.i!=c.j) spent[c.j]+=c.to_j;
        skip[index]=1;
        ++out.statistics.skipped_quartets;
    }
    out.statistics.max_error_bound=static_cast<double>(*std::max_element(spent.begin(),spent.end()));
    const auto t2=Clock::now();
    out.statistics.selection_seconds=std::chrono::duration<double>(t2-t1).count();

    std::vector<LocalPreparedShellPair> prepared;
    prepared.reserve(np);
    for (const auto& p:pairs)
        prepared.push_back(prepare_local_shell_pair(shells[p.a],shells[p.b]));
    const auto t3=Clock::now();
    out.statistics.pair_preparation_seconds=std::chrono::duration<double>(t3-t2).count();
    Eigen::MatrixXd upper=Eigen::MatrixXd::Zero(nao,nao);
    std::size_t index=0;
    for (std::size_t ip=0;ip<np;++ip) {
        const auto [a,b]=std::pair{pairs[ip].a,pairs[ip].b};
        for (std::size_t iq=ip;iq<np;++iq,++index) {
            if (skip[index]) continue;
            const auto [c,d]=std::pair{pairs[iq].a,pairs[iq].b};
            const auto q=local_shell_quartet_paired(prepared[ip],prepared[iq]);
            detail::accumulate_local_density_quartet(upper, density, offsets, a, b, c, d, ip == iq, q);
        }
    }
    out.matrix=Eigen::MatrixXd(upper.selfadjointView<Eigen::Upper>());
    out.statistics.contraction_seconds=std::chrono::duration<double>(Clock::now()-t3).count();
    return out;
}

} // namespace chem::kumo
