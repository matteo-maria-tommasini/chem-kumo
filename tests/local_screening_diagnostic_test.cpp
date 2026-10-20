#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
struct Pair {
    std::size_t a, b;
    double bmax, weighted_density;
    long double work;
};
struct Candidate {
    std::size_t i, j;
    double to_i, to_j;
    long double work;
};
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        using namespace chem::kumo;
        for (int arg = 1; arg < argc; ++arg) {
            const auto molecule = Molecule::from_fchk(argv[arg]);
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const Eigen::MatrixXd norms = local_ao_product_norms(basis);
            std::vector<Eigen::Index> offsets(shells.size()+1, 0);
            std::vector<long double> primitive_counts(shells.size(), 0.0L);
            for (std::size_t a = 0; a < shells.size(); ++a) {
                offsets[a+1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
                for (const auto& component : shells[a])
                    primitive_counts[a] += static_cast<long double>(component.size());
            }
            for (Spin spin : {Spin::Total, Spin::Difference}) {
                const Eigen::MatrixXd density = molecule.density_matrix(spin);
                if (density.rows() != norms.rows() || density.cols() != norms.cols() || !density.allFinite())
                    throw std::runtime_error("density shape or finiteness mismatch");
                const double asym = (density-density.transpose()).cwiseAbs().maxCoeff();
                if (asym > 1e-10*std::max(1.0,density.cwiseAbs().maxCoeff()))
                    throw std::runtime_error("asymmetric density");
                std::vector<Pair> pairs;
                for (std::size_t a=0; a<shells.size(); ++a)
                    for (std::size_t b=a; b<shells.size(); ++b) {
                        double bmax=0.0, dsum=0.0;
                        for (Eigen::Index mu=offsets[a]; mu<offsets[a+1]; ++mu)
                            for (Eigen::Index nu=offsets[b]; nu<offsets[b+1]; ++nu) {
                                const double v=norms(mu,nu);
                                bmax=std::max(bmax,v);
                                dsum+=std::abs(density(mu,nu))*v;
                            }
                        pairs.push_back({a,b,bmax,dsum,primitive_counts[a]*primitive_counts[b]});
                    }
                std::vector<Candidate> candidates;
                long double total_work=0.0L;
                for (std::size_t i=0; i<pairs.size(); ++i)
                    for (std::size_t j=i; j<pairs.size(); ++j) {
                        const auto& x=pairs[i];
                        const auto& y=pairs[j];
                        const double mi=x.a==x.b ? 1.0 : 2.0;
                        const double mj=y.a==y.b ? 1.0 : 2.0;
                        const double ei=mj*x.bmax*y.weighted_density;
                        const double ej=(i==j) ? 0.0 : mi*y.bmax*x.weighted_density;
                        const long double work=x.work*y.work;
                        candidates.push_back({i,j,ei,ej,work});
                        total_work+=work;
                    }
                // Cheap-to-discard quartets first. Deterministic tie-break by pair indices.
                std::sort(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b){
                    const double aa=std::max(a.to_i,a.to_j);
                    const double bb=std::max(b.to_i,b.to_j);
                    if (aa!=bb) return aa<bb;
                    if (a.i!=b.i) return a.i<b.i;
                    return a.j<b.j;
                });
                std::cout<<std::scientific<<std::setprecision(12)
                         <<"system="<<argv[arg]
                         <<" spin="<<(spin==Spin::Total?"total":"difference")
                         <<" nAO="<<norms.rows()<<" shell_pairs="<<pairs.size()
                         <<" quartets="<<candidates.size()<<'\n';
                for (double epsilon : {1e-6,1e-8,1e-10,1e-12}) {
                    std::vector<long double> spent(pairs.size(),0.0L);
                    std::size_t skipped=0;
                    long double saved=0.0L;
                    for (const auto& c : candidates) {
                        if (spent[c.i]+static_cast<long double>(c.to_i)>epsilon) continue;
                        if (c.i!=c.j && spent[c.j]+static_cast<long double>(c.to_j)>epsilon) continue;
                        spent[c.i]+=c.to_i;
                        if (c.i!=c.j) spent[c.j]+=c.to_j;
                        ++skipped;
                        saved+=c.work;
                    }
                    const long double max_bound=*std::max_element(spent.begin(),spent.end());
                    if (max_bound>epsilon*(1.0L+1e-12L))
                        throw std::runtime_error("screening budget exceeded");
                    std::cout<<" epsilon="<<epsilon
                             <<" skipped_quartets="<<skipped
                             <<" skipped_fraction="<<static_cast<double>(skipped)/candidates.size()
                             <<" max_error_bound="<<static_cast<double>(max_bound)
                             <<" estimated_saved_work_fraction="<<static_cast<double>(saved/total_work)
                             <<'\n';
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr<<"screening diagnostic: "<<e.what()<<'\n';
        return 1;
    }
}
