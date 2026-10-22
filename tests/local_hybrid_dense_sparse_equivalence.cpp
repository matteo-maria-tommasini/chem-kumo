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
            std::vector<Eigen::MatrixXd> cs_screened, hybrid_screened, cs_spent, hybrid_spent;
            double bound_seconds=0, cs_selection_seconds=0, hybrid_selection_seconds=0;
            double cs_contraction_seconds=0, hybrid_contraction_seconds=0;
            std::size_t cs_evaluated=0, hybrid_evaluated=0;
            for (std::size_t k=0;k<densities.size();++k) {
                cs_screened.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                hybrid_screened.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                cs_spent.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                hybrid_spent.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
            }
            detail::LocalSparseBudget cs_local(nao,budgets), hybrid_local(nao,budgets);
            double cs_sparse_seconds=0, hybrid_sparse_seconds=0;
            std::size_t decision_mismatches=0, quartets=0;
            std::size_t cs_skipped=0, hybrid_skipped=0;
            for (std::size_t x=0;x<pairs.size();++x)
                for (std::size_t y=x;y<pairs.size();++y) {
                    const auto& a=pairs[x];
                    const auto& b=pairs[y];
                    ++quartets;
                    std::vector<Eigen::MatrixXd> cs_charge, hybrid_charge;
                    std::vector<detail::LocalBudgetCharge> cs_local_charges, hybrid_local_charges;
                    for (std::size_t k=0;k<densities.size();++k) {
                        cs_charge.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                        hybrid_charge.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                    }
                    const auto bound_start=Clock::now();
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
                                    const double spatial=local_shell_pair_envelope_bound(
                                        a.env,i,j,b.env,l,m);
                                    if (!(cs>=0) || !(spatial>=0) || !std::isfinite(cs) || !std::isfinite(spatial))
                                        throw std::runtime_error("nonfinite or negative bound");
                                    const long double bound=std::min(cs,spatial);

                                    for (std::size_t k=0;k<densities.size();++k) {
                                        const long double direct=mb*std::abs(static_cast<long double>(densities[k](lambda,sigma)));
                                        const double dc=static_cast<double>(direct*cs);
                                        const double dh=static_cast<double>(direct*bound);
                                        cs_charge[k](mu,nu)+=dc;
                                        hybrid_charge[k](mu,nu)+=dh;
                                        if (dc!=0.0) cs_local_charges.push_back({mu,nu,k,dc});
                                        if (dh!=0.0) hybrid_local_charges.push_back({mu,nu,k,dh});
                                        if (x!=y && lambda<=sigma) {
                                            const long double factor=(a.a==a.b && mu!=nu)?2.0L:ma;
                                            const long double reciprocal=factor*std::abs(static_cast<long double>(densities[k](mu,nu)));
                                            const double rc=static_cast<double>(reciprocal*cs);
                                            const double rh=static_cast<double>(reciprocal*bound);
                                            cs_charge[k](lambda,sigma)+=rc;
                                            hybrid_charge[k](lambda,sigma)+=rh;
                                            if (rc!=0.0) cs_local_charges.push_back({lambda,sigma,k,rc});
                                            if (rh!=0.0) hybrid_local_charges.push_back({lambda,sigma,k,rh});
                                        }
                                    }
                                }
                        }
                    auto admissible = [&](const std::vector<Eigen::MatrixXd>& spent,
                                          const std::vector<Eigen::MatrixXd>& charge) {
                        for (std::size_t k=0;k<densities.size();++k)
                            for (Eigen::Index i=0;i<nao;++i)
                                for (Eigen::Index j=i;j<nao;++j)
                                    if (static_cast<long double>(spent[k](i,j)) + charge[k](i,j) > budgets[k])
                                        return false;
                        return true;
                    };
                    bound_seconds+=elapsed(bound_start,Clock::now());
                    const auto cs_sel_start=Clock::now();
                    const bool omit_cs=admissible(cs_spent,cs_charge);
                    cs_selection_seconds+=elapsed(cs_sel_start,Clock::now());
                    const auto hybrid_sel_start=Clock::now();
                    const bool omit_hybrid=admissible(hybrid_spent,hybrid_charge);
                    hybrid_selection_seconds+=elapsed(hybrid_sel_start,Clock::now());
                    const auto cs_sparse_start=Clock::now();
                    const bool sparse_cs=cs_local.try_charge(cs_local_charges);
                    cs_sparse_seconds+=elapsed(cs_sparse_start,Clock::now());
                    const auto hybrid_sparse_start=Clock::now();
                    const bool sparse_hybrid=hybrid_local.try_charge(hybrid_local_charges);
                    hybrid_sparse_seconds+=elapsed(hybrid_sparse_start,Clock::now());
                    if (omit_cs!=sparse_cs || omit_hybrid!=sparse_hybrid) {
                        ++decision_mismatches;
                        throw std::runtime_error("dense/sparse decision mismatch at quartet " + std::to_string(quartets));
                    }
                    if (omit_cs) {
                        ++cs_skipped;
                        for (std::size_t k=0;k<densities.size();++k) cs_spent[k]+=cs_charge[k];
                    }
                    if (omit_hybrid) {
                        ++hybrid_skipped;
                        for (std::size_t k=0;k<densities.size();++k) hybrid_spent[k]+=hybrid_charge[k];
                    }
                    // Only evaluate a quartet if at least one strategy retains it.
                    if (!omit_cs || !omit_hybrid) {
                        const auto q=local_shell_quartet({&shells[a.a],&shells[a.b],&shells[b.a],&shells[b.b]});
                        // A shared evaluation is charged to both strategies conservatively.
                        // Separate timed contractions exclude the common integral construction.
                        if (!omit_cs) {
                            ++cs_evaluated;
                            const auto t=Clock::now();
                            for (std::size_t k=0;k<densities.size();++k)
                                detail::accumulate_local_density_quartet(cs_screened[k],densities[k],offsets,
                                    a.a,a.b,b.a,b.b,x==y,q);
                            cs_contraction_seconds+=elapsed(t,Clock::now());
                        }
                        if (!omit_hybrid) {
                            ++hybrid_evaluated;
                            const auto t=Clock::now();
                            for (std::size_t k=0;k<densities.size();++k)
                                detail::accumulate_local_density_quartet(hybrid_screened[k],densities[k],offsets,
                                    a.a,a.b,b.a,b.b,x==y,q);
                            hybrid_contraction_seconds+=elapsed(t,Clock::now());
                        }
                    }
                }
            double max_budget_difference=0.0;
            for (std::size_t k=0;k<densities.size();++k) {
                max_budget_difference=std::max(max_budget_difference,
                    (cs_spent[k]-cs_local.spent()[k]).cwiseAbs().maxCoeff());
                max_budget_difference=std::max(max_budget_difference,
                    (hybrid_spent[k]-hybrid_local.spent()[k]).cwiseAbs().maxCoeff());
            }
            if (max_budget_difference > 1e-18)
                throw std::runtime_error("dense/sparse budget ledger mismatch");
            // This is a selection/contraction diagnostic, not an end-to-end speedup:
            // the common integral evaluations and the shared bound construction are
            // intentionally not attributed to either strategy.
            std::cout<<std::scientific<<std::setprecision(12)
                     <<"system="<<argv[arg]
                     <<" shell_pairs="<<pairs.size()<<" quartets="<<quartets
                     <<" cs_skipped="<<cs_skipped<<" hybrid_skipped="<<hybrid_skipped
                     <<" cs_evaluated="<<cs_evaluated<<" hybrid_evaluated="<<hybrid_evaluated
                     <<" norms_seconds="<<norms_seconds
                     <<" envelopes_seconds="<<envelope_seconds
                     <<" shared_bounds_seconds="<<bound_seconds
                     <<" cs_selection_seconds="<<cs_selection_seconds
                     <<" hybrid_selection_seconds="<<hybrid_selection_seconds
                     <<" cs_sparse_seconds="<<cs_sparse_seconds
                     <<" hybrid_sparse_seconds="<<hybrid_sparse_seconds
                     <<" cs_selection_speedup="<<(cs_sparse_seconds>0?cs_selection_seconds/cs_sparse_seconds:0)
                     <<" hybrid_selection_speedup="<<(hybrid_sparse_seconds>0?hybrid_selection_seconds/hybrid_sparse_seconds:0)
                     <<" decision_mismatches="<<decision_mismatches
                     <<" max_budget_difference="<<max_budget_difference
                     <<" cs_accumulation_seconds="<<cs_contraction_seconds
                     <<" hybrid_accumulation_seconds="<<hybrid_contraction_seconds
                     <<" status=OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
