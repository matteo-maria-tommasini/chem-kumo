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
#include <limits>
#include <iostream>
#include <numeric>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using Clock=std::chrono::steady_clock;
    auto elapsed=[](auto a, auto b){return std::chrono::duration<double>(b-a).count();};
    using namespace chem::kumo;
    int repeats=9;
    int first=1;
    if (argc >= 3 && std::string(argv[1]) == "--repeats") {
        repeats=std::stoi(argv[2]);
        first=3;
    }
    if (repeats<3 || repeats>101 || argc<=first) {
        std::cerr << "usage: benchmark [--repeats 9] file.fchk[.gz] ...\n";
        return 2;
    }
    auto median=[](std::vector<double> v) {
        std::sort(v.begin(),v.end());
        return v[v.size()/2];
    };
    try {
        for (int arg=first; arg<argc; ++arg) {
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
            // Warm up both paths; alternate execution order to reduce cache/order bias.
            (void)run(false);
            (void)run(true);
            std::vector<double> cs_compute_times, hybrid_compute_times;
            std::vector<double> cs_integral_times, hybrid_integral_times;
            std::vector<double> cs_selection_times, hybrid_selection_times;
            std::size_t cs_skipped=0, hybrid_skipped=0;
            double worst_difference=0, checksum=0;
            for (int rep=0;rep<repeats;++rep) {
                Result cs, hy;
                if (rep%2==0) { cs=run(false); hy=run(true); }
                else { hy=run(true); cs=run(false); }
                if (rep==0) { cs_skipped=cs.skipped; hybrid_skipped=hy.skipped; }
                if (cs.skipped!=cs_skipped || hy.skipped!=hybrid_skipped ||
                    hy.skipped<cs.skipped)
                    throw std::runtime_error("non-reproducible screening decisions");
                for (std::size_t k=0;k<densities.size();++k)
                    worst_difference=std::max(worst_difference,
                        (cs.w[k]-hy.w[k]).cwiseAbs().maxCoeff());
                if (!std::isfinite(cs.checksum) || !std::isfinite(hy.checksum))
                    throw std::runtime_error("nonfinite checksum");
                checksum+=cs.checksum+hy.checksum;
                cs_compute_times.push_back(cs.selection+cs.integrals+cs.accumulation);
                hybrid_compute_times.push_back(hy.selection+hy.integrals+hy.accumulation);
                cs_integral_times.push_back(cs.integrals);
                hybrid_integral_times.push_back(hy.integrals);
                cs_selection_times.push_back(cs.selection);
                hybrid_selection_times.push_back(hy.selection);
            }
            const double tc=median(cs_compute_times), th=median(hybrid_compute_times);
            const double extra_setup=envelope_seconds; // Shared bounds contain BOTH methods: cannot attribute their cost.
            const double saved=tc-th;
            const double break_even=saved>0?extra_setup/saved:std::numeric_limits<double>::infinity();
            std::cout<<std::scientific<<std::setprecision(12)
                <<"system="<<argv[arg]
                <<" repeats="<<repeats
                <<" quartets="<<candidates.size()
                <<" cs_skipped="<<cs_skipped<<" hybrid_skipped="<<hybrid_skipped
                <<" cs_integrals_median="<<median(cs_integral_times)
                <<" hybrid_integrals_median="<<median(hybrid_integral_times)
                <<" cs_selection_median="<<median(cs_selection_times)
                <<" hybrid_selection_median="<<median(hybrid_selection_times)
                <<" cs_compute_median="<<tc
                <<" hybrid_compute_median="<<th
                <<" compute_speedup="<<(th>0?tc/th:0)
                <<" extra_envelope_setup_seconds="<<extra_setup
                <<" shared_bounds_seconds="<<shared_bounds_seconds
                <<" norms_seconds="<<norms_seconds
                <<" break_even_reuses_lower_bound="<<break_even
                <<" max_cs_hybrid_difference="<<worst_difference
                <<" checksum="<<checksum
                <<" status=OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
