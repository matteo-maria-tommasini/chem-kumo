#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/local_density_quartet_contraction.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/local_sparse_budget.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using Clock=std::chrono::steady_clock;
    auto elapsed=[](auto a, auto b){return std::chrono::duration<double>(b-a).count();};
    using namespace chem::kumo;
    if (argc < 2) { std::cerr << "usage: test file.fchk[.gz] ...\n"; return 2; }
    try {
        for (int arg=1; arg<argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule=Molecule::from_fchk(argv[arg]);
            const auto shells=local_shell_components(basis);
            std::vector<Eigen::Index> offsets(shells.size()+1,0);
            for (std::size_t s=0;s<shells.size();++s)
                offsets[s+1]=offsets[s]+static_cast<Eigen::Index>(shells[s].size());
            const Eigen::Index nao=offsets.back();
            const std::vector<Eigen::MatrixXd> densities{
                molecule.density_matrix(Spin::Total), molecule.density_matrix(Spin::Difference)};
            const std::vector<double> budgets{1e-8,1e-10};
            for (const auto& d : densities)
                if (d.rows()!=nao || d.cols()!=nao || !d.allFinite())
                    throw std::runtime_error("invalid physical density");
            const auto norm_start=Clock::now();
            const auto norms=local_ao_product_norms(basis);
            const double norms_seconds=elapsed(norm_start,Clock::now());
            const std::size_t ns=std::min<std::size_t>(4,shells.size());
            std::vector<std::size_t> chosen;
            for (std::size_t t=0;t<ns;++t) {
                const auto s=ns==1?0:t*(shells.size()-1)/(ns-1);
                if (chosen.empty() || chosen.back()!=s) chosen.push_back(s);
            }
            struct Pair { std::size_t a,b; LocalShellPairEnvelope env; };
            const auto envelope_start=Clock::now();
            std::vector<Pair> pairs;
            for (std::size_t a=0;a<chosen.size();++a)
                for (std::size_t b=a;b<chosen.size();++b)
                    pairs.push_back({chosen[a],chosen[b],
                        prepare_local_shell_pair_envelope(shells[chosen[a]],shells[chosen[b]],0.7)});
            const double envelope_seconds=elapsed(envelope_start,Clock::now());
            struct Candidate {
                std::size_t x, y;
                std::vector<detail::LocalBudgetCharge> cs, hybrid;
            };
            std::vector<Candidate> candidates;
            const auto bound_start=Clock::now();
            for (std::size_t x=0;x<pairs.size();++x)
                for (std::size_t y=x;y<pairs.size();++y) {
                    const auto& a=pairs[x];
                    const auto& b=pairs[y];
                    Candidate candidate{x,y,{},{}};
                    const long double mb=b.a==b.b?1.0L:2.0L;
                    const long double ma=a.a==a.b?1.0L:2.0L;
                    for (std::size_t i=0;i<a.env.na;++i)
                        for (std::size_t j=0;j<a.env.nb;++j) {
                            const auto mu=offsets[a.a]+static_cast<Eigen::Index>(i);
                            const auto nu=offsets[a.b]+static_cast<Eigen::Index>(j);
                            if (mu>nu) continue;
                            for (std::size_t l=0;l<b.env.na;++l)
                                for (std::size_t m=0;m<b.env.nb;++m) {
                                    const auto lambda=offsets[b.a]+static_cast<Eigen::Index>(l);
                                    const auto sigma=offsets[b.b]+static_cast<Eigen::Index>(m);
                                    const double cs=norms(mu,nu)*norms(lambda,sigma);
                                    const double spatial=local_shell_pair_envelope_bound(a.env,i,j,b.env,l,m);
                                    if (!std::isfinite(cs) || !std::isfinite(spatial) || cs<0 || spatial<0)
                                        throw std::runtime_error("invalid bound");
                                    const double hybrid=std::min(cs,spatial);
                                    for (std::size_t k=0;k<densities.size();++k) {
                                        const long double direct=mb*std::abs(static_cast<long double>(densities[k](lambda,sigma)));
                                        const double dc=static_cast<double>(direct*cs);
                                        const double dh=static_cast<double>(direct*hybrid);
                                        if (dc!=0) candidate.cs.push_back({mu,nu,k,dc});
                                        if (dh!=0) candidate.hybrid.push_back({mu,nu,k,dh});
                                        if (x!=y && lambda<=sigma) {
                                            const long double factor=(a.a==a.b && mu!=nu)?2.0L:ma;
                                            const long double reciprocal=factor*std::abs(static_cast<long double>(densities[k](mu,nu)));
                                            const double rc=static_cast<double>(reciprocal*cs);
                                            const double rh=static_cast<double>(reciprocal*hybrid);
                                            if (rc!=0) candidate.cs.push_back({lambda,sigma,k,rc});
                                            if (rh!=0) candidate.hybrid.push_back({lambda,sigma,k,rh});
                                        }
                                    }
                                }
                        }
                    candidates.push_back(std::move(candidate));
                }
            const double shared_bounds_seconds=elapsed(bound_start,Clock::now());
            struct Result {
                std::vector<Eigen::MatrixXd> w;
                std::size_t skipped=0, evaluated=0;
                double selection=0, integrals=0, accumulation=0;
                double checksum=0;
            };
            auto run=[&](bool hybrid) {
                Result result;
                for (std::size_t k=0;k<densities.size();++k)
                    result.w.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                detail::LocalSparseBudget ledger(nao,budgets);
                for (const auto& candidate:candidates) {
                    const auto ts=Clock::now();
                    const bool skip=ledger.try_charge(hybrid?candidate.hybrid:candidate.cs);
                    result.selection+=elapsed(ts,Clock::now());
                    if (skip) { ++result.skipped; continue; }
                    const auto& a=pairs[candidate.x];
                    const auto& b=pairs[candidate.y];
                    const auto ti=Clock::now();
                    const auto q=local_shell_quartet({&shells[a.a],&shells[a.b],&shells[b.a],&shells[b.b]});
                    result.integrals+=elapsed(ti,Clock::now());
                    const auto ta=Clock::now();
                    for (std::size_t k=0;k<densities.size();++k)
                        detail::accumulate_local_density_quartet(result.w[k],densities[k],offsets,
                            a.a,a.b,b.a,b.b,candidate.x==candidate.y,q);
                    result.accumulation+=elapsed(ta,Clock::now());
                    ++result.evaluated;
                }
                for (const auto& w:result.w) result.checksum+=w.squaredNorm();
                return result;
            };
            // Independent contractions: integrals are generated only for retained quartets.
            const auto cs=run(false);
            const auto hybrid=run(true);
            if (hybrid.skipped<cs.skipped)
                throw std::runtime_error("hybrid skipped fewer quartets");
            double max_difference=0;
            for (std::size_t k=0;k<densities.size();++k)
                max_difference=std::max(max_difference,(cs.w[k]-hybrid.w[k]).cwiseAbs().maxCoeff());
            const double cs_compute=cs.selection+cs.integrals+cs.accumulation;
            const double hybrid_compute=hybrid.selection+hybrid.integrals+hybrid.accumulation;
            const double cs_amortized=norms_seconds+cs_compute;
            const double hybrid_amortized=norms_seconds+envelope_seconds+shared_bounds_seconds+hybrid_compute;
            std::cout<<std::scientific<<std::setprecision(12)
                <<"system="<<argv[arg]
                <<" quartets="<<candidates.size()
                <<" cs_skipped="<<cs.skipped<<" hybrid_skipped="<<hybrid.skipped
                <<" cs_evaluated="<<cs.evaluated<<" hybrid_evaluated="<<hybrid.evaluated
                <<" norms_seconds="<<norms_seconds
                <<" envelopes_seconds="<<envelope_seconds
                <<" shared_bounds_seconds="<<shared_bounds_seconds
                <<" cs_selection_seconds="<<cs.selection
                <<" hybrid_selection_seconds="<<hybrid.selection
                <<" cs_integrals_seconds="<<cs.integrals
                <<" hybrid_integrals_seconds="<<hybrid.integrals
                <<" cs_accumulation_seconds="<<cs.accumulation
                <<" hybrid_accumulation_seconds="<<hybrid.accumulation
                <<" compute_speedup="<<(hybrid_compute>0?cs_compute/hybrid_compute:0)
                <<" amortized_speedup="<<(hybrid_amortized>0?cs_amortized/hybrid_amortized:0)
                <<" max_cs_hybrid_difference="<<max_difference
                <<" cs_checksum="<<cs.checksum<<" hybrid_checksum="<<hybrid.checksum
                <<" status=OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
