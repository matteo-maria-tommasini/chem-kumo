#pragma once

#include "chem/kumo/local_density_quartet_contraction.hpp"

#include "chem/kumo/local_density_overlap_screened.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

enum class ScreeningOrder { Canonical, Importance };

struct LocalStreamingStatistics {
    std::size_t shell_pairs = 0;
    std::size_t total_quartets = 0;
    std::size_t skipped_quartets = 0;
    std::size_t prepared_pair_instances = 0;
    std::size_t peak_pair_bytes = 0;
    double max_error_bound = 0.0;
    double norm_seconds = 0.0;
    double ordering_seconds = 0.0;
    double contraction_seconds = 0.0;
};

struct LocalStreamingDensityOverlap {
    Eigen::MatrixXd matrix;
    LocalStreamingStatistics statistics;
};

// Two active prepared shell pairs (outer reused, inner prepared on demand).
// memory_limit_bytes caps the measured pair-vector storage, not the total RSS.
// No O(N_shell_pairs^2) candidate list, no full Gaussian-pair cache.
[[nodiscard]] inline LocalStreamingDensityOverlap local_density_overlap_streaming(
    const BasisSet& basis, const Eigen::Ref<const Eigen::MatrixXd>& density,
    double epsilon, std::size_t memory_limit_bytes,
    ScreeningOrder order = ScreeningOrder::Canonical) {
    using Clock = std::chrono::steady_clock;
    if (!(epsilon >= 0.0) || !std::isfinite(epsilon))
        throw std::invalid_argument("streaming overlap: invalid epsilon");
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size()+1,0);
    for (std::size_t a=0;a<shells.size();++a)
        offsets[a+1]=offsets[a]+static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao=offsets.back();
    if (density.rows()!=nao || density.cols()!=nao || !density.allFinite())
        throw std::invalid_argument("streaming overlap: invalid density");
    LocalStreamingDensityOverlap out;
    out.matrix=Eigen::MatrixXd::Zero(nao,nao);
    if (!nao) return out;
    const double scale=density.cwiseAbs().maxCoeff();
    if ((density-density.transpose()).cwiseAbs().maxCoeff() >
        1e-10*std::max(1.0,scale))
        throw std::invalid_argument("streaming overlap: asymmetric density");
    struct Pair { std::size_t a,b; long double bmax=0, dsum=0; };
    std::vector<Pair> pairs;
    for (std::size_t a=0;a<shells.size();++a)
        for (std::size_t b=a;b<shells.size();++b)
            pairs.push_back({a,b});
    const std::size_t np=pairs.size();
    out.statistics.shell_pairs=np;
    out.statistics.total_quartets=np*(np+1)/2;
    if (scale==0.0) {
        out.statistics.skipped_quartets=out.statistics.total_quartets;
        return out;
    }
    const auto start=Clock::now();
    const Eigen::MatrixXd norms=local_ao_product_norms(basis);
    out.statistics.norm_seconds=std::chrono::duration<double>(Clock::now()-start).count();
    for (auto& p:pairs)
        for (Eigen::Index mu=offsets[p.a];mu<offsets[p.a+1];++mu)
            for (Eigen::Index nu=offsets[p.b];nu<offsets[p.b+1];++nu) {
                const long double b=norms(mu,nu);
                p.bmax=std::max(p.bmax,b);
                p.dsum+=std::abs(static_cast<long double>(density(mu,nu)))*b;
            }
    // Sort indices, never the prepared Gaussian intermediates.
    const auto ordering_start=Clock::now();
    std::vector<std::size_t> permutation(np);
    for (std::size_t i=0;i<np;++i) permutation[i]=i;
    if (order==ScreeningOrder::Importance) {
        const auto importance=[&](std::size_t i) {
            const auto& p=pairs[i];
            return (p.a==p.b?1.0L:2.0L)*p.bmax*p.dsum;
        };
        std::sort(permutation.begin(),permutation.end(),[&](std::size_t i,std::size_t j) {
            const auto si=importance(i), sj=importance(j);
            return si<sj || (si==sj && i<j);
        });
    }
    out.statistics.ordering_seconds=std::chrono::duration<double>(Clock::now()-ordering_start).count();
    std::vector<long double> spent(np,0.0L);
    Eigen::MatrixXd upper=Eigen::MatrixXd::Zero(nao,nao);
    // Bound the vector-owned allocations, including their capacities and headers.
    const auto bytes=[](const LocalPreparedShellPair& p)->std::size_t {
        std::size_t n=p.components.capacity()*sizeof(std::vector<GaussianPairIntermediate>);
        for (const auto& c:p.components) {
            const auto cap=c.capacity();
            if (cap>(std::numeric_limits<std::size_t>::max()-n)/sizeof(GaussianPairIntermediate))
                throw std::overflow_error("streaming overlap: byte count overflow");
            n+=cap*sizeof(GaussianPairIntermediate);
        }
        return n;
    };
    const auto estimate_pair_bytes=[&](const Pair& p)->std::size_t {
        const auto& a=shells[p.a];
        const auto& b=shells[p.b];
        std::size_t count=0;
        for (const auto& ca:a)
            for (const auto& cb:b) {
                if (ca.size() && cb.size()> (std::numeric_limits<std::size_t>::max()-count)/ca.size())
                    throw std::overflow_error("streaming overlap: primitive count overflow");
                count+=ca.size()*cb.size();
            }
        const std::size_t headers=a.size()*b.size()*sizeof(std::vector<GaussianPairIntermediate>);
        if (count>(std::numeric_limits<std::size_t>::max()-headers)/sizeof(GaussianPairIntermediate))
            throw std::overflow_error("streaming overlap: pair estimate overflow");
        return headers+count*sizeof(GaussianPairIntermediate);
    };
    const auto contraction_start=Clock::now();
    for (std::size_t pi=0;pi<np;++pi) {
        const std::size_t ip=permutation[pi];
        const auto& x=pairs[ip];
        const long double mi=x.a==x.b?1.0L:2.0L;
        // Delay preparing the outer pair until the first surviving quartet.
        LocalPreparedShellPair outer;
        bool outer_ready=false;
        for (std::size_t pj=pi;pj<np;++pj) {
            const std::size_t iq=permutation[pj];
            const auto& y=pairs[iq];
            const long double mj=y.a==y.b?1.0L:2.0L;
            const long double to_i=mj*x.bmax*y.dsum;
            const long double to_j=ip==iq?0.0L:mi*y.bmax*x.dsum;
            if (spent[ip]+to_i<=static_cast<long double>(epsilon) &&
                (ip==iq || spent[iq]+to_j<=static_cast<long double>(epsilon))) {
                spent[ip]+=to_i;
                if (ip!=iq) spent[iq]+=to_j;
                ++out.statistics.skipped_quartets;
                continue;
            }
            const std::size_t estimated=estimate_pair_bytes(x)
                +(ip==iq?0:estimate_pair_bytes(y));
            if (estimated>memory_limit_bytes)
                throw std::runtime_error("streaming overlap: Gaussian-pair memory limit exceeded");
            if (!outer_ready) {
                outer=prepare_local_shell_pair(shells[x.a],shells[x.b]);
                ++out.statistics.prepared_pair_instances;
                outer_ready=true;
            }
            LocalPreparedShellPair inner;
            const LocalPreparedShellPair* second=&outer;
            if (ip!=iq) {
                inner=prepare_local_shell_pair(shells[y.a],shells[y.b]);
                ++out.statistics.prepared_pair_instances;
                second=&inner;
            }
            const std::size_t resident=bytes(outer)+(ip==iq?0:bytes(inner));
            out.statistics.peak_pair_bytes=std::max(out.statistics.peak_pair_bytes,resident);
            if (resident>memory_limit_bytes)
                throw std::runtime_error("streaming overlap: Gaussian-pair memory limit exceeded");
            const auto q=local_shell_quartet_paired(outer,*second);
            detail::accumulate_local_density_quartet(upper, density, offsets, x.a, x.b, y.a, y.b, ip == iq, q);
        }
    }
    out.matrix=Eigen::MatrixXd(upper.selfadjointView<Eigen::Upper>());
    out.statistics.max_error_bound=static_cast<double>(*std::max_element(spent.begin(),spent.end()));
    out.statistics.contraction_seconds=std::chrono::duration<double>(Clock::now()-contraction_start).count();
    return out;
}

} // namespace chem::kumo
