#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/local_quartet_sampling.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/local_sparse_budget.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using namespace chem::kumo;
using Charge = detail::LocalBudgetCharge;
using Key = std::tuple<std::size_t, Eigen::Index, Eigen::Index>;
using Ledger = std::map<Key, long double>;

// Build one transaction per strategy, including both canonical contraction directions.
// Long-double accumulation reduces cancellation/rounding loss; final conversion rounds up.
void add(Ledger& ledger, std::size_t channel, Eigen::Index row, Eigen::Index col,
         long double amount) {
    if (row > col) throw std::runtime_error("noncanonical AO charge");
    if (!(amount >= 0) || !std::isfinite(amount)) throw std::runtime_error("invalid AO charge");
    ledger[{channel,row,col}] += amount;
    if (!std::isfinite(ledger[{channel,row,col}])) throw std::runtime_error("charge overflow");
}
std::vector<Charge> materialize(const Ledger& ledger) {
    std::vector<Charge> result;
    result.reserve(ledger.size());
    for (const auto& [key, amount] : ledger) {
        double rounded = static_cast<double>(amount);
        if (!std::isfinite(rounded)) throw std::runtime_error("charge conversion overflow");
        if (static_cast<long double>(rounded) < amount)
            rounded = std::nextafter(rounded, std::numeric_limits<double>::infinity());
        if (!std::isfinite(rounded)) throw std::runtime_error("charge upward rounding overflow");
        result.push_back({std::get<1>(key),std::get<2>(key),std::get<0>(key),rounded});
    }
    return result;
}
}

int main(int argc, char** argv) {
    using namespace chem::kumo;
    using Clock = std::chrono::steady_clock;
    try {
        std::size_t target = 1000;
        int arg = 1;
        if (argc > 2 && std::string(argv[1]) == "--target") {
            target = std::stoull(argv[2]); arg = 3;
        }
        if (arg >= argc || target == 0 || target > 10000 || target % 4)
            throw std::invalid_argument("usage: profile [--target 1000|5000|10000] file.fchk[.gz] ...");
        std::cout << std::scientific << std::setprecision(10);
        for (; arg < argc; ++arg) {
            const std::string path = argv[arg];
            const BasisSet basis{FormattedCheckpoint(path)};
            const auto molecule = Molecule::from_fchk(path);
            const auto shells = local_shell_components(basis);
            const auto t0 = Clock::now();
            const auto sampled = sampling::sample_quartets(shells, target / 4);
            const auto t1 = Clock::now();
            const auto norms = local_ao_product_norms(basis);
            const auto t2 = Clock::now();
            std::vector<Eigen::Index> offsets(shells.size()+1,0);
            for (std::size_t s=0;s<shells.size();++s)
                offsets[s+1]=offsets[s]+static_cast<Eigen::Index>(shells[s].size());
            const auto nao=offsets.back();
            const std::vector<Eigen::MatrixXd> densities{
                molecule.density_matrix(Spin::Total),
                molecule.density_matrix(Spin::Difference)};
            for (const auto& d : densities)
                if (d.rows()!=nao || d.cols()!=nao || !d.allFinite())
                    throw std::runtime_error("invalid physical density");
            const std::vector<double> eps{1e-8,1e-10};
            detail::LocalSparseBudget cs_budget(nao,eps), hybrid_budget(nao,eps);
            std::map<std::size_t,LocalShellPairEnvelope> cache;
            auto envelope = [&](std::size_t p)->const LocalShellPairEnvelope& {
                auto it=cache.find(p);
                if(it==cache.end()) {
                    const auto& sp=sampled.pairs.at(p);
                    it=cache.emplace(p,prepare_local_shell_pair_envelope(
                        shells[sp.a],shells[sp.b],0.7)).first;
                }
                return it->second;
            };
            std::array<std::size_t,4> cs_skip{},hybrid_skip{},invalid{},ao_elements{};
            std::array<long double,4> total_work{},cs_work{},hybrid_work{};
            const auto t3=Clock::now();
            for (unsigned bin=0;bin<4;++bin) {
                const auto& samples=sampled.samples[bin];
                for (const auto& sample : samples) {
                    const auto [pi,pj]=sampling::decode(sample.id);
                    const auto& ab=sampled.pairs.at(pi);
                    const auto& cd=sampled.pairs.at(pj);
                    const auto& ea=envelope(pi);
                    const auto& eb=envelope(pj);
                    Ledger cs_ledger, hybrid_ledger;
                    const long double mult_cd=(cd.a==cd.b)?1.0L:2.0L;
                    const long double mult_ab=(ab.a==ab.b)?1.0L:2.0L;
                    for (std::size_t i=0;i<ea.na;++i)
                        for (std::size_t j=0;j<ea.nb;++j) {
                            const auto mu=offsets[ab.a]+static_cast<Eigen::Index>(i);
                            const auto nu=offsets[ab.b]+static_cast<Eigen::Index>(j);
                            if(mu>nu) continue;
                            for (std::size_t l=0;l<eb.na;++l)
                                for (std::size_t m=0;m<eb.nb;++m) {
                                    const auto lambda=offsets[cd.a]+static_cast<Eigen::Index>(l);
                                    const auto sigma=offsets[cd.b]+static_cast<Eigen::Index>(m);
                                    const double cs=norms(mu,nu)*norms(lambda,sigma);
                                    const double spatial=local_shell_pair_envelope_bound(ea,i,j,eb,l,m);
                                    if (!(cs>=0) || !(spatial>=0) || !std::isfinite(cs) || !std::isfinite(spatial))
                                        throw std::runtime_error("invalid bound");
                                    const long double hybrid=std::min(cs,spatial);
                                    ++ao_elements[bin];
                                    for (std::size_t k=0;k<densities.size();++k) {
                                        const long double direct=mult_cd*std::abs(static_cast<long double>(densities[k](lambda,sigma)));
                                        add(cs_ledger,k,mu,nu,direct*cs);
                                        add(hybrid_ledger,k,mu,nu,direct*hybrid);
                                        if(pi!=pj && lambda<=sigma) {
                                            const long double factor=(ab.a==ab.b && mu!=nu)?2.0L:mult_ab;
                                            const long double reciprocal=factor*std::abs(static_cast<long double>(densities[k](mu,nu)));
                                            add(cs_ledger,k,lambda,sigma,reciprocal*cs);
                                            add(hybrid_ledger,k,lambda,sigma,reciprocal*hybrid);
                                        }
                                    }
                                }
                        }
                    const auto cs_charges=materialize(cs_ledger);
                    const auto hybrid_charges=materialize(hybrid_ledger);
                    // Strictly zero budgets cannot authorize nonempty quartets: conservative policy.
                    const bool allow=std::all_of(eps.begin(),eps.end(),[](double e){return e>0;});
                    const bool omit_cs=allow && cs_budget.try_charge(cs_charges);
                    const bool omit_hybrid=allow && hybrid_budget.try_charge(hybrid_charges);
                    total_work[bin]+=sample.cost;
                    if(omit_cs) {++cs_skip[bin];cs_work[bin]+=sample.cost;}
                    if(omit_hybrid) {++hybrid_skip[bin];hybrid_work[bin]+=sample.cost;}
                }
                std::cout<<"system="<<path<<" target="<<target<<" bin="<<bin
                         <<" sampled="<<samples.size()<<" ao_elements="<<ao_elements[bin]
                         <<" invalid="<<invalid[bin]<<" cs_skipped="<<cs_skip[bin]
                         <<" hybrid_skipped="<<hybrid_skip[bin]
                         <<" work_proxy="<<total_work[bin]
                         <<" cs_saved_work_proxy="<<cs_work[bin]
                         <<" hybrid_saved_work_proxy="<<hybrid_work[bin]<<'\n';
            }
            const auto t4=Clock::now();
            for (std::size_t k=0;k<eps.size();++k) {
                const double cs_max=cs_budget.spent()[k].maxCoeff();
                const double hybrid_max=hybrid_budget.spent()[k].maxCoeff();
                if(cs_max>eps[k] || hybrid_max>eps[k]) throw std::runtime_error("budget exceeded");
                std::cout<<"system="<<path<<" channel="<<k<<" epsilon="<<eps[k]
                         <<" cs_max_budget_used="<<cs_max
                         <<" hybrid_max_budget_used="<<hybrid_max<<'\n';
            }
            std::cout<<"system="<<path
                     <<" sample_seconds="<<std::chrono::duration<double>(t1-t0).count()
                     <<" norms_seconds="<<std::chrono::duration<double>(t2-t1).count()
                     <<" screening_seconds="<<std::chrono::duration<double>(t4-t3).count()
                     <<" envelope_cache_entries="<<cache.size()
                     <<" mode=density_weighted_full_AO_sample status=OK\n";
        }
    } catch(const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
