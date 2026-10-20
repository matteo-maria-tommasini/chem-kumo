#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Pair {
    std::size_t a, b;
    long double bmax = 0.0L;
    long double dsum = 0.0L;
    long double work = 0.0L;
};
struct Result {
    std::size_t skipped = 0;
    long double max_bound = 0.0L;
    long double saved_work = 0.0L;
    long double total_work = 0.0L;
    double seconds = 0.0;
};

// Reference streaming selector (same canonical traversal as Stage 03C.11d).
// Deterministic, strictly streaming: O(number of shell pairs) auxiliary memory.
Result select_streaming(const std::vector<Pair>& pairs, long double epsilon) {
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    const std::size_t n = pairs.size();
    std::vector<long double> spent(n, 0.0L);
    Result out;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& x = pairs[i];
        const long double mi = x.a == x.b ? 1.0L : 2.0L;
        for (std::size_t j = i; j < n; ++j) {
            const auto& y = pairs[j];
            const long double mj = y.a == y.b ? 1.0L : 2.0L;
            const long double to_i = mj * x.bmax * y.dsum;
            const long double to_j = i == j ? 0.0L : mi * y.bmax * x.dsum;
            const long double work = x.work * y.work;
            out.total_work += work;
            if (spent[i] + to_i > epsilon ||
                (i != j && spent[j] + to_j > epsilon)) continue;
            spent[i] += to_i;
            if (i != j) spent[j] += to_j;
            ++out.skipped;
            out.saved_work += work;
        }
    }
    if (!spent.empty()) out.max_bound = *std::max_element(spent.begin(), spent.end());
    out.seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return out;
}

// Stage 03C.12b: segment tree with deferred normalized error budgets.
// r_i=(epsilon-E_i)/b_i, so a reciprocal update subtracts a constant
// from a contiguous range of r_i values. b_i=0 uses +infinity.
// The traversal preserves the exact canonical (i,j) greedy order.
struct LazyTree {
    std::size_t n;
    std::vector<long double> minimum, lazy, asum, wsum;
    explicit LazyTree(const std::vector<Pair>& p, long double eps)
        : n(p.size()), minimum(4*n+4), lazy(4*n+4),
          asum(4*n+4), wsum(4*n+4) {
        if (n) build(1,0,n,p,eps);
    }
    void build(std::size_t node,std::size_t l,std::size_t r,
               const std::vector<Pair>& p,long double eps) {
        if(r-l==1) {
            const auto& q=p[l];
            minimum[node]=q.bmax>0 ? eps/q.bmax
                : std::numeric_limits<long double>::infinity();
            asum[node]=(q.a==q.b ? 1.0L:2.0L)*q.dsum;
            wsum[node]=q.work;
            return;
        }
        auto m=l+(r-l)/2;
        build(node*2,l,m,p,eps); build(node*2+1,m,r,p,eps);
        pull(node);
    }
    void pull(std::size_t k) {
        minimum[k]=std::min(minimum[k*2],minimum[k*2+1]);
        asum[k]=asum[k*2]+asum[k*2+1];
        wsum[k]=wsum[k*2]+wsum[k*2+1];
    }
    void apply(std::size_t k,long double v) {
        if(std::isfinite(minimum[k])) minimum[k]-=v;
        lazy[k]+=v;
    }
    void push(std::size_t k) {
        if(lazy[k]!=0) {
            apply(k*2,lazy[k]); apply(k*2+1,lazy[k]);
            lazy[k]=0;
        }
    }
    void subtract(std::size_t k,std::size_t l,std::size_t r,
                  std::size_t ql,std::size_t qr,long double v) {
        if(qr<=l||r<=ql||v==0) return;
        if(ql<=l&&r<=qr) { apply(k,v); return; }
        push(k); auto m=l+(r-l)/2;
        subtract(k*2,l,m,ql,qr,v);
        subtract(k*2+1,m,r,ql,qr,v);
        pull(k);
    }
    void subtract(std::size_t l,std::size_t r,long double v) {
        if(l<r) subtract(1,0,n,l,r,v);
    }
    long double point(std::size_t k,std::size_t l,std::size_t r,std::size_t i) {
        if(r-l==1) return minimum[k];
        push(k); auto m=l+(r-l)/2;
        return i<m ? point(k*2,l,m,i) : point(k*2+1,m,r,i);
    }
    long double point(std::size_t i) {return point(1,0,n,i);}
    long double min_all() const {return minimum[1];}
    long double sum_a(std::size_t k,std::size_t l,std::size_t r,
                      std::size_t ql,std::size_t qr) const {
        if(qr<=l||r<=ql) return 0;
        if(ql<=l&&r<=qr) return asum[k];
        auto m=l+(r-l)/2;
        return sum_a(k*2,l,m,ql,qr)+sum_a(k*2+1,m,r,ql,qr);
    }
    long double sum_a(std::size_t l,std::size_t r) const {
        return sum_a(1,0,n,l,r);
    }
    long double sum_work(std::size_t k,std::size_t l,std::size_t r,
                         std::size_t ql,std::size_t qr) const {
        if(qr<=l||r<=ql) return 0;
        if(ql<=l&&r<=qr) return wsum[k];
        auto m=l+(r-l)/2;
        return sum_work(k*2,l,m,ql,qr)+sum_work(k*2+1,m,r,ql,qr);
    }
    long double sum_work(std::size_t l,std::size_t r) const {
        return sum_work(1,0,n,l,r);
    }
};
struct LazyResult : Result {
    std::size_t aggregate_decisions=0, individual_visits=0;
    std::size_t accepted_by_aggregate=0, node_visits=0;
    std::size_t tree_bytes=0;
};
LazyResult select_lazy(const std::vector<Pair>& pairs,long double eps) {
    using Clock=std::chrono::steady_clock;
    const auto start=Clock::now();
    LazyResult out;
    const std::size_t n=pairs.size();
    if(!n) return out;
    LazyTree tree(pairs,eps);
    out.tree_bytes=(tree.minimum.capacity()+tree.lazy.capacity()+
                    tree.asum.capacity()+tree.wsum.capacity())*sizeof(long double);
    const long double total_pair_work=tree.sum_work(0,n);
    out.total_work=(total_pair_work*total_pair_work+
        [&] {long double v=0; for(const auto& p:pairs) v+=p.work*p.work; return v;}())/2;
    for(std::size_t i=0;i<n;++i) {
        const long double ai=(pairs[i].a==pairs[i].b?1.0L:2.0L)*pairs[i].dsum;

        // Handle diagonal separately: only one update, not two.
        ++out.individual_visits;
        if(tree.point(i)>=ai) {
            tree.subtract(i,i+1,ai);
            ++out.skipped;
            out.saved_work+=pairs[i].work*pairs[i].work;
        }
        if(i+1==n) continue;
        // Visit j>i in canonical order; aggregate accepted nodes are
        // updated lazily, while rejected nodes are recursively subdivided.
        auto visit=[&](auto&& self,std::size_t k,std::size_t l,std::size_t r)->void {
            if(r<=i+1) return;
            ++out.node_visits;
            if(l>i && tree.point(i)>=tree.asum[k] && tree.minimum[k]>=ai) {
                // The two inequalities guarantee every leaf budget.
                tree.subtract(l,r,ai);
                tree.subtract(i,i+1,tree.asum[k]);
                out.skipped+=r-l;
                out.saved_work+=pairs[i].work*tree.wsum[k];
                ++out.aggregate_decisions;
                out.accepted_by_aggregate+=r-l;
                return;
            }
            if(r-l==1) {
                ++out.individual_visits;
                if(tree.point(i)>=tree.asum[k] && tree.minimum[k]>=ai) {
                    tree.subtract(l,r,ai);
                    tree.subtract(i,i+1,tree.asum[k]);
                    ++out.skipped;
                    out.saved_work+=pairs[i].work*pairs[l].work;
                }
                return;
            }
            tree.push(k);
            auto m=l+(r-l)/2;
            self(self,k*2,l,m);
            self(self,k*2+1,m,r);
            tree.pull(k);
        };
        visit(visit,1,0,n);
    }
    // Reconstruct E_i for the diagnostic; no per-quartet state is retained.
    long double max_bound=0;
    for(std::size_t i=0;i<n;++i) if(pairs[i].bmax>0) {
        const long double spent=eps-pairs[i].bmax*tree.point(i);
        max_bound=std::max(max_bound,spent);
    }
    out.max_bound=max_bound;
    out.seconds=std::chrono::duration<double>(Clock::now()-start).count();
    return out;
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kumo-local-screening-lazy-test <fchk> [<fchk> ...]\n";
        return 2;
    }
    try {
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = 1; arg < argc; ++arg) {
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const Eigen::MatrixXd norms = local_ao_product_norms(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            std::vector<long double> primitives(shells.size(), 0.0L);
            for (std::size_t a = 0; a < shells.size(); ++a) {
                offsets[a+1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
                for (const auto& component : shells[a]) primitives[a] += component.size();
            }
            for (Spin spin : {Spin::Total, Spin::Difference}) {
                const Eigen::MatrixXd density = molecule.density_matrix(spin);
                if (density.rows() != norms.rows() || density.cols() != norms.cols() ||
                    !density.allFinite()) throw std::runtime_error("invalid density matrix");
                if ((density-density.transpose()).cwiseAbs().maxCoeff() >
                    1e-10 * std::max(1.0, density.cwiseAbs().maxCoeff()))
                    throw std::runtime_error("density matrix not symmetric");
                std::vector<Pair> pairs;
                pairs.reserve(shells.size() * (shells.size() + 1) / 2);
                for (std::size_t a = 0; a < shells.size(); ++a)
                    for (std::size_t b = a; b < shells.size(); ++b) {
                        Pair p{a, b};
                        p.work = primitives[a] * primitives[b];
                        for (Eigen::Index mu = offsets[a]; mu < offsets[a+1]; ++mu)
                            for (Eigen::Index nu = offsets[b]; nu < offsets[b+1]; ++nu) {
                                const long double v = norms(mu, nu);
                                p.bmax = std::max(p.bmax, v);
                                p.dsum += std::abs(static_cast<long double>(density(mu,nu))) * v;
                            }
                        pairs.push_back(p);
                    }
                const std::size_t n = pairs.size();
                std::cout << "system=" << argv[arg]
                          << " spin=" << (spin == Spin::Total ? "total" : "difference")
                          << " nAO=" << norms.rows() << " shell_pairs=" << n
                          << " quartets=" << n*(n+1)/2
                          << " reference_aux_bytes=" << n*sizeof(long double) << '\n';
                for (long double epsilon : {1e-6L, 1e-8L, 1e-10L, 1e-12L}) {
                    const Result ref = select_streaming(pairs, epsilon);
                    const LazyResult r=select_lazy(pairs,epsilon);
                    const long double fraction=r.total_work>0 ? r.saved_work/r.total_work:0;
                    if(r.skipped!=ref.skipped ||
                       std::abs(r.max_bound-ref.max_bound)>1e-11L*std::max(1.0L,epsilon))
                        throw std::runtime_error("lazy/reference mismatch");
                    if(r.max_bound>epsilon+1e-17L)
                        throw std::runtime_error("lazy budget exceeded");
                    std::cout << " epsilon=" << static_cast<double>(epsilon)
                              << " reference_skipped=" << ref.skipped
                              << " lazy_skipped=" << r.skipped
                              << " aggregate_decisions=" << r.aggregate_decisions
                              << " aggregate_quartets=" << r.accepted_by_aggregate
                              << " individual_visits=" << r.individual_visits
                              << " node_visits=" << r.node_visits
                              << " tree_bytes=" << r.tree_bytes
                              << " max_error_bound=" << static_cast<double>(r.max_bound)
                              << " saved_work_fraction=" << static_cast<double>(fraction)
                              << " selection_seconds=" << r.seconds << '\n';
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "lazy screening diagnostic: " << e.what() << '\n';
        return 1;
    }
}
