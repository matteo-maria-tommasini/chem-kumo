#pragma once

#include "chem/kumo/local_density_quartet_contraction.hpp"

#include "chem/kumo/local_density_overlap_screened.hpp"
#include <chrono>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace chem::kumo {

/// Immutable geometry-dependent preparation for repeated density contractions.
/// A plan is tied to the BasisSet geometry used at construction time.
/// Do not mutate the basis or share a plan with a different AO ordering.
class LocalDensityOverlapPlan {
public:
    explicit LocalDensityOverlapPlan(const BasisSet& basis) {
        using Clock = std::chrono::steady_clock;
        const auto t0 = Clock::now();
        const auto shells = local_shell_components(basis);
        offsets_.assign(shells.size()+1, 0);
        for (std::size_t a=0; a<shells.size(); ++a)
            offsets_[a+1] = offsets_[a] + static_cast<Eigen::Index>(shells[a].size());
        for (std::size_t a=0; a<shells.size(); ++a)
            for (std::size_t b=a; b<shells.size(); ++b)
                pairs_.push_back({a,b,0.0});
        const auto t1 = Clock::now();
        norms_ = local_ao_product_norms(basis);
        const auto t2 = Clock::now();
        for (auto& pair : pairs_)
            for (Eigen::Index mu=offsets_[pair.a]; mu<offsets_[pair.a+1]; ++mu)
                for (Eigen::Index nu=offsets_[pair.b]; nu<offsets_[pair.b+1]; ++nu)
                    pair.bmax=std::max(pair.bmax,norms_(mu,nu));
        prepared_.reserve(pairs_.size());
        for (const auto& p : pairs_)
            prepared_.push_back(prepare_local_shell_pair(shells[p.a],shells[p.b]));
        const auto t3 = Clock::now();
        norm_seconds_ = std::chrono::duration<double>(t2-t1).count();
        preparation_seconds_ = std::chrono::duration<double>(t3-t2).count()
            + std::chrono::duration<double>(t1-t0).count();
    }

    [[nodiscard]] Eigen::Index ao_count() const noexcept { return offsets_.back(); }
    [[nodiscard]] std::size_t shell_pair_count() const noexcept { return pairs_.size(); }
    [[nodiscard]] double norm_preparation_seconds() const noexcept { return norm_seconds_; }
    [[nodiscard]] double pair_preparation_seconds() const noexcept { return preparation_seconds_; }
    [[nodiscard]] std::size_t primitive_pair_count() const {
        std::size_t n=0;
        for (const auto& p:prepared_) n+=p.pair_count();
        return n;
    }

    [[nodiscard]] LocalScreenedDensityOverlap evaluate(
        const Eigen::Ref<const Eigen::MatrixXd>& density, double epsilon) const {
        using Clock = std::chrono::steady_clock;
        if (!std::isfinite(epsilon) || epsilon < 0.0)
            throw std::invalid_argument("LocalDensityOverlapPlan: invalid error budget");
        const Eigen::Index nao=ao_count();
        if (density.rows()!=nao || density.cols()!=nao || !density.allFinite())
            throw std::invalid_argument("LocalDensityOverlapPlan: invalid density");
        LocalScreenedDensityOverlap out;
        out.matrix=Eigen::MatrixXd::Zero(nao,nao);
        const std::size_t np=pairs_.size();
        out.statistics.total_quartets=np*(np+1)/2;
        if (!nao) return out;
        const double scale=density.cwiseAbs().maxCoeff();
        if ((density-density.transpose()).cwiseAbs().maxCoeff() >
            1e-10*std::max(1.0,scale))
            throw std::invalid_argument("LocalDensityOverlapPlan: asymmetric density");
        if (scale==0.0) {
            out.statistics.skipped_quartets=out.statistics.total_quartets;
            return out;
        }
        const auto t0=Clock::now();
        std::vector<long double> dsum(np,0.0L);
        for (std::size_t ip=0; ip<np; ++ip) {
            const auto& pair=pairs_[ip];
            for (Eigen::Index mu=offsets_[pair.a];mu<offsets_[pair.a+1];++mu)
                for (Eigen::Index nu=offsets_[pair.b];nu<offsets_[pair.b+1];++nu)
                    dsum[ip]+=static_cast<long double>(std::abs(density(mu,nu)))
                             *static_cast<long double>(norms_(mu,nu));
        }
        struct Candidate { std::size_t i,j; long double to_i,to_j; };
        std::vector<Candidate> candidates;
        candidates.reserve(out.statistics.total_quartets);
        for (std::size_t i=0;i<np;++i)
            for (std::size_t j=i;j<np;++j) {
                const auto& x=pairs_[i];
                const auto& y=pairs_[j];
                const long double mi=x.a==x.b?1.0L:2.0L;
                const long double mj=y.a==y.b?1.0L:2.0L;
                candidates.push_back({i,j,
                    mj*static_cast<long double>(x.bmax)*dsum[j],
                    i==j?0.0L:mi*static_cast<long double>(y.bmax)*dsum[i]});
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
        out.statistics.max_error_bound=static_cast<double>(
            *std::max_element(spent.begin(),spent.end()));
        const auto t1=Clock::now();
        out.statistics.selection_seconds=std::chrono::duration<double>(t1-t0).count();
        Eigen::MatrixXd upper=Eigen::MatrixXd::Zero(nao,nao);
        std::size_t index=0;
        for (std::size_t ip=0;ip<np;++ip) {
            const auto [a,b]=std::pair{pairs_[ip].a,pairs_[ip].b};
            for (std::size_t iq=ip;iq<np;++iq,++index) {
                if (skip[index]) continue;
                const auto [c,d]=std::pair{pairs_[iq].a,pairs_[iq].b};
                const auto q=local_shell_quartet_paired(prepared_[ip],prepared_[iq]);
                detail::accumulate_local_density_quartet(upper, density, offsets_, a, b, c, d, ip == iq, q);
            }
        }
        out.matrix=Eigen::MatrixXd(upper.selfadjointView<Eigen::Upper>());
        out.statistics.contraction_seconds=
            std::chrono::duration<double>(Clock::now()-t1).count();
        return out;
    }

private:
    struct Pair { std::size_t a,b; double bmax; };
    std::vector<Eigen::Index> offsets_;
    std::vector<Pair> pairs_;
    Eigen::MatrixXd norms_;
    std::vector<LocalPreparedShellPair> prepared_;
    double norm_seconds_=0.0;
    double preparation_seconds_=0.0;
};

} // namespace chem::kumo
