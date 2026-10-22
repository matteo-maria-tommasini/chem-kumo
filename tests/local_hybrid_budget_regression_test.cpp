#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/local_density_quartet_contraction.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
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
            const auto norms=local_ao_product_norms(basis);
            const std::size_t ns=std::min<std::size_t>(4,shells.size());
            std::vector<std::size_t> chosen;
            for (std::size_t t=0;t<ns;++t) {
                const auto s=ns==1?0:t*(shells.size()-1)/(ns-1);
                if (chosen.empty() || chosen.back()!=s) chosen.push_back(s);
            }
            struct Pair { std::size_t a,b; LocalShellPairEnvelope env; };
            std::vector<Pair> pairs;
            for (std::size_t a=0;a<chosen.size();++a)
                for (std::size_t b=a;b<chosen.size();++b)
                    pairs.push_back({chosen[a],chosen[b],
                        prepare_local_shell_pair_envelope(shells[chosen[a]],shells[chosen[b]],0.7)});
            std::vector<Eigen::MatrixXd> exact, screened, spent;
            for (std::size_t k=0;k<densities.size();++k) {
                exact.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                screened.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
                spent.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
            }
            std::size_t skipped=0, retained=0, hybrid_better=0, elements=0;
            for (std::size_t x=0;x<pairs.size();++x)
                for (std::size_t y=x;y<pairs.size();++y) {
                    const auto& a=pairs[x];
                    const auto& b=pairs[y];
                    std::vector<Eigen::MatrixXd> charge;
                    for (std::size_t k=0;k<densities.size();++k)
                        charge.emplace_back(Eigen::MatrixXd::Zero(nao,nao));
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
                                    hybrid_better+=(spatial<cs);
                                    ++elements;
                                    for (std::size_t k=0;k<densities.size();++k) {
                                        charge[k](mu,nu)+=static_cast<double>(mb*std::abs(
                                            static_cast<long double>(densities[k](lambda,sigma)))*bound);
                                        if (x!=y && lambda<=sigma) {
                                            const long double factor=(a.a==a.b && mu!=nu)?2.0L:ma;
                                            charge[k](lambda,sigma)+=static_cast<double>(factor*std::abs(
                                                static_cast<long double>(densities[k](mu,nu)))*bound);
                                        }
                                    }
                                }
                        }
                    bool omit=true;
                    for (std::size_t k=0;k<densities.size() && omit;++k)
                        for (Eigen::Index i=0;i<nao && omit;++i)
                            for (Eigen::Index j=i;j<nao;++j)
                                if (static_cast<long double>(spent[k](i,j)) + charge[k](i,j) > budgets[k]) {
                                    omit=false;
                                    break;
                                }
                    if (omit) {
                        ++skipped;
                        for (std::size_t k=0;k<densities.size();++k) spent[k]+=charge[k];
                    } else ++retained;
                    const auto q=local_shell_quartet({&shells[a.a],&shells[a.b],&shells[b.a],&shells[b.b]});
                    for (std::size_t k=0;k<densities.size();++k) {
                        detail::accumulate_local_density_quartet(exact[k],densities[k],offsets,
                            a.a,a.b,b.a,b.b,x==y,q);
                        if (!omit)
                            detail::accumulate_local_density_quartet(screened[k],densities[k],offsets,
                                a.a,a.b,b.a,b.b,x==y,q);
                    }
                }
            std::cout<<std::scientific<<std::setprecision(12);
            for (std::size_t k=0;k<densities.size();++k) {
                exact[k]=Eigen::MatrixXd(exact[k].selfadjointView<Eigen::Upper>());
                screened[k]=Eigen::MatrixXd(screened[k].selfadjointView<Eigen::Upper>());
                const double err=(exact[k]-screened[k]).cwiseAbs().maxCoeff();
                const double bound=spent[k].cwiseAbs().maxCoeff();
                const double scale=std::max(1.0,exact[k].cwiseAbs().maxCoeff());
                if (err>bound+2e-11*scale || bound>budgets[k]+1e-18)
                    throw std::runtime_error("hybrid cumulative error budget violated");
                std::cout<<"system="<<argv[arg]<<" channel="<<k
                         <<" max_abs="<<err<<" max_bound="<<bound
                         <<" epsilon="<<budgets[k]<<'\n';
            }
            std::cout<<"system="<<argv[arg]<<" shell_pairs="<<pairs.size()
                     <<" quartets="<<skipped+retained<<" skipped="<<skipped
                     <<" retained="<<retained<<" hybrid_tighter="<<hybrid_better
                     <<" ao_elements="<<elements<<" status=OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
