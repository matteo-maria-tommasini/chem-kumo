#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
struct Pair { long double bmax=0, work=0; std::vector<long double> d; };
struct Fenwick {
    std::vector<long double> counts, weights;
    explicit Fenwick(std::size_t n): counts(n+1,0), weights(n+1,0) {}
    void add(std::size_t pos, long double weight) {
        for (++pos; pos<counts.size(); pos+=pos&-pos) { counts[pos]+=1; weights[pos]+=weight; }
    }
    std::pair<long double,long double> prefix(std::size_t count) const {
        long double c=0,w=0;
        for (;count;count-=count&-count) {c+=counts[count];w+=weights[count];}
        return {c,w};
    }
};

// A directed suffix is independently budget-certified for every output pair.
// Its mutual intersection is a subset of every directed suffix, hence also certified.
struct MutualResult {
    std::size_t skipped=0, nodes=0, queries=0;
    long double saved_work=0, total_work=0, max_budget_ratio=0;
    double seconds=0;
};
MutualResult mutual_tails(const std::vector<Pair>& p,
                          const std::vector<long double>& eps,
                          const std::vector<std::size_t>& order,
                          bool verify) {
    const auto start=Clock::now();
    const std::size_t n=p.size(), k=eps.size();
    std::vector<std::vector<long double>> suffix(k,std::vector<long double>(n+1,0));
    for(std::size_t t=n;t-->0;)
        for(std::size_t c=0;c<k;++c)
            suffix[c][t]=suffix[c][t+1]+p[order[t]].d[c];
    std::vector<std::size_t> cut(n,n);
    MutualResult out;
    for(std::size_t i=0;i<n;++i) {
        std::size_t lo=0,hi=n;
        while(lo<hi) {
            const auto mid=lo+(hi-lo)/2;
            bool fits=true;
            for(std::size_t c=0;c<k;++c) {
                if(p[order[i]].bmax*suffix[c][mid]>eps[c]*(1.0L-1e-12L)) {
                    fits=false;break;
                }
            }
            ++out.nodes;
            if(fits) hi=mid; else lo=mid+1;
        }
        cut[i]=lo;
        ++out.queries;
    }
    // Process rank i ascending. Active ranks j have cut[j]<=i.
    // Query active j>=cut[i], then remove diagonal and divide by two.
    std::vector<std::vector<std::size_t>> activate(n+1);
    for(std::size_t j=0;j<n;++j) activate[cut[j]].push_back(j);
    Fenwick tree(n);
    long double ordered_count=0, ordered_weight=0, diagonal_count=0, diagonal_weight=0;
    long double total_pair_work=0, sum_work=0, sum_work_squared=0;
    for(const auto& x:p) {sum_work+=x.work;sum_work_squared+=x.work*x.work;}
    total_pair_work=(sum_work*sum_work+sum_work_squared)/2;
    for(std::size_t i=0;i<n;++i) {
        for(auto j:activate[i]) if(j<n) tree.add(j,p[order[j]].work);
        const auto [all_c,all_w]=tree.prefix(n);
        const auto [before_c,before_w]=tree.prefix(cut[i]);
        ordered_count+=all_c-before_c;
        ordered_weight+=p[order[i]].work*(all_w-before_w);
        if(i>=cut[i]) {++diagonal_count;diagonal_weight+=p[order[i]].work*p[order[i]].work;}
    }
    out.skipped=static_cast<std::size_t>((ordered_count+diagonal_count)/2);
    out.saved_work=(ordered_weight+diagonal_weight)/2;
    out.total_work=total_pair_work;
    // Explicit verification only for small systems, including actual reciprocal budgets.
    if(verify && n<=1024) {
        std::size_t count=0;
        long double weight=0;
        std::vector<std::vector<long double>> spent(k,std::vector<long double>(n,0));
        for(std::size_t i=0;i<n;++i) for(std::size_t j=i;j<n;++j) {
            if(j<cut[i] || i<cut[j]) continue;
            ++count;weight+=p[order[i]].work*p[order[j]].work;
            for(std::size_t c=0;c<k;++c) {
                spent[c][i]+=p[order[i]].bmax*p[order[j]].d[c];
                if(i!=j) spent[c][j]+=p[order[j]].bmax*p[order[i]].d[c];
            }
        }
        if(count!=out.skipped || std::abs(weight-out.saved_work)>1e-8L*std::max(1.0L,weight))
            throw std::runtime_error("mutual intersection count/weight mismatch");
        for(std::size_t c=0;c<k;++c) for(std::size_t i=0;i<n;++i) {
            if(spent[c][i]>eps[c]+1e-18L*std::max(1.0L,eps[c]))
                throw std::runtime_error("mutual budget exceeded");
            if(eps[c]>0) out.max_budget_ratio=std::max(out.max_budget_ratio,spent[c][i]/eps[c]);
        }
    }
    out.seconds=std::chrono::duration<double>(Clock::now()-start).count();
    return out;
}
}
int main(int argc,char**argv) {
    if(argc<2) {std::cerr<<"usage: diagnostic <fchk> [<fchk> ...]\n";return 2;}
    try {
        using namespace chem::kumo;
        std::cout<<std::scientific<<std::setprecision(12);
        for(int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto mol=Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd total=mol.density_matrix(Spin::Total);
            const Eigen::MatrixXd spin=mol.density_matrix(Spin::Difference);
            const std::vector<Eigen::MatrixXd> densities{0.5*(total+spin),0.5*(total-spin)};
            const std::vector<long double> budgets{1e-8L,1e-10L};
            const auto shells=local_shell_components(basis);
            const auto norm_start=Clock::now();
            const auto norms=local_ao_product_norms(basis);
            const double norm_seconds=std::chrono::duration<double>(Clock::now()-norm_start).count();
            std::vector<Eigen::Index> offsets(shells.size()+1,0);
            std::vector<long double> primitive(shells.size(),0);
            for(std::size_t a=0;a<shells.size();++a) {
                offsets[a+1]=offsets[a]+static_cast<Eigen::Index>(shells[a].size());
                for(const auto& component:shells[a]) primitive[a]+=component.size();
            }
            if(offsets.back()!=total.rows()) throw std::runtime_error("AO count mismatch");
            std::vector<Pair> pairs;
            for(std::size_t a=0;a<shells.size();++a) for(std::size_t b=a;b<shells.size();++b) {
                Pair p; p.work=primitive[a]*primitive[b]; p.d.assign(2,0);
                const long double mult=(a==b?1.0L:2.0L);
                for(Eigen::Index mu=offsets[a];mu<offsets[a+1];++mu)
                    for(Eigen::Index nu=offsets[b];nu<offsets[b+1];++nu) {
                        const long double v=norms(mu,nu);
                        p.bmax=std::max(p.bmax,v);
                        for(std::size_t c=0;c<2;++c)
                            p.d[c]+=mult*std::abs(static_cast<long double>(densities[c](mu,nu)))*v;
                    }
                pairs.push_back(std::move(p));
            }
            const auto n=pairs.size();
            std::vector<std::size_t> order(n);std::iota(order.begin(),order.end(),0);
            std::vector<long double> scores(n,0);
            for(std::size_t i=0;i<n;++i) for(std::size_t c=0;c<2;++c)
                scores[i]=std::max(scores[i],pairs[i].bmax*pairs[i].d[c]/budgets[c]);
            std::sort(order.begin(),order.end(),[&](auto a,auto b){
                return scores[a]>scores[b] || (scores[a]==scores[b] && a<b);
            });
            const auto result=mutual_tails(pairs,budgets,order,true);
            const long double formal=static_cast<long double>(n)*(n+1)/2;
            std::cout<<"system="<<argv[arg]<<" nao="<<total.rows()<<" shell_pairs="<<n
                     <<" formal_quartets="<<std::fixed<<std::setprecision(0)<<formal
                     <<std::scientific<<std::setprecision(12)
                     <<" norm_seconds="<<norm_seconds
                     <<" queries="<<result.queries<<" binary_nodes="<<result.nodes
                     <<" mutual_skipped="<<result.skipped
                     <<" mutual_skip_fraction="<<static_cast<long double>(result.skipped)/formal
                     <<" estimated_saved_work_fraction="<<(result.total_work>0?result.saved_work/result.total_work:0)
                     <<" max_budget_ratio="<<result.max_budget_ratio<<" selection_seconds="<<result.seconds<<"\n";
        }
    } catch(const std::exception& e) {std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
