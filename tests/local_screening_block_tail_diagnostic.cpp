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
struct Result { std::size_t directed_skips=0, queries=0, nodes=0; long double saved_work=0, total_work=0; double seconds=0; };
// For every fixed output shell pair i, omit a suffix of *input* shell pairs.
// This is a one-way diagnostic, NOT a replacement for symmetric quartet screening.
Result select_tails(const std::vector<Pair>& p, const std::vector<long double>& eps,
                    const std::vector<std::size_t>& order, bool verify) {
    const std::size_t n=p.size(), k=eps.size();
    std::vector<std::vector<long double>> suffix(k, std::vector<long double>(n+1,0));
    std::vector<long double> work_suffix(n+1,0);
    for (std::size_t t=n; t-->0;) {
        const auto& x=p[order[t]];
        work_suffix[t]=work_suffix[t+1]+x.work;
        for(std::size_t c=0;c<k;++c) suffix[c][t]=suffix[c][t+1]+x.d[c];
    }
    const auto start=Clock::now();
    Result out;
    for (std::size_t i=0;i<n;++i) {
        const long double b=p[i].bmax;
        // Smallest index t such that all channels satisfy b*sum_{j>=t} d_j <= eps.
        std::size_t lo=0, hi=n;
        while(lo<hi) {
            const auto mid=lo+(hi-lo)/2;
            bool fits=true;
            for(std::size_t c=0;c<k;++c) {
                // Conservative numerical margin on the allowed budget.
                const long double safe=eps[c]*(1.0L-1e-12L);
                if(b*suffix[c][mid]>safe) { fits=false; break; }
            }
            ++out.nodes;
            if(fits) hi=mid; else lo=mid+1;
        }
        const auto cut=lo;
        if(verify && n<=1024) {
            for(std::size_t c=0;c<k;++c) {
                long double direct=0;
                for(std::size_t t=cut;t<n;++t) direct+=p[order[t]].d[c];
                if(b*direct>eps[c]) throw std::runtime_error("invalid directed tail bound");
            }
        }
        out.directed_skips+=n-cut;
        out.saved_work+=p[i].work*work_suffix[cut];
        out.total_work+=p[i].work*work_suffix[0];
        ++out.queries;
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
            const auto result=select_tails(pairs,budgets,order,true);
            const long double formal=static_cast<long double>(n)*n;
            std::cout<<"system="<<argv[arg]<<" nao="<<total.rows()<<" shell_pairs="<<n
                     <<" formal_directed="<<std::fixed<<std::setprecision(0)<<formal
                     <<std::scientific<<std::setprecision(12)
                     <<" norm_seconds="<<norm_seconds
                     <<" queries="<<result.queries<<" binary_nodes="<<result.nodes
                     <<" directed_tail_skips="<<result.directed_skips
                     <<" directed_skip_fraction="<<static_cast<long double>(result.directed_skips)/formal
                     <<" estimated_saved_work_fraction="<<(result.total_work>0?result.saved_work/result.total_work:0)
                     <<" selection_seconds="<<result.seconds<<"\n";
        }
    } catch(const std::exception& e) {std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
