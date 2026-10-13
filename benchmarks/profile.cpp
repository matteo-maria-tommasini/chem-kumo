#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
using Clock = std::chrono::steady_clock;
static double seconds(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b-a).count();
}
static Eigen::MatrixXd unpack(const std::vector<double>& v, Eigen::Index n) {
    if (v.size()!=static_cast<std::size_t>(n*(n+1)/2)) throw std::runtime_error("bad density matrix");
    Eigen::MatrixXd p(n,n);
    std::size_t k=0;
    for(Eigen::Index i=0;i<n;++i) for(Eigen::Index j=0;j<=i;++j) {
        p(i,j)=p(j,i)=v[k++];
    }
    return p;
}
int main(int argc, char** argv) {
    if(argc<2 || argc>5) {
        std::cerr<<"usage: kumo-profile file.fchk [spacing=0.3] [repeats=3] [AO_threshold=1e-8]\n";
        return 2;
    }
    try {
        const double spacing=argc>2?std::stod(argv[2]):0.3;
        const int repeats=argc>3?std::stoi(argv[3]):3;
        const double threshold=argc>4?std::stod(argv[4]):1e-8;
        if(!(spacing>0 && std::isfinite(spacing)) || repeats<1 || repeats>50 ||
           !(threshold>=0 && std::isfinite(threshold))) throw std::invalid_argument("invalid arguments");
        const FormattedCheckpoint f(argv[1]);
        const BasisEvaluator evaluator{BasisSet(f)};
        const auto mol=Molecule::from_fchk(argv[1]);
        const Grid grid=Grid::around(mol.atoms(),3.0,spacing);
        const Eigen::Index n=static_cast<Eigen::Index>(evaluator.basis_size());
        const auto p=unpack(f.reals("Total SCF Density"),n);
        const std::size_t count=static_cast<std::size_t>(grid.counts[0])*grid.counts[1]*grid.counts[2];
        std::vector<Eigen::VectorXd> aos;
        aos.reserve(count);
        double ao_best=1e100;
        for(int rep=0;rep<repeats;++rep) {
            aos.clear();
            const auto t0=Clock::now();
            for(int i=0;i<grid.counts[0];++i)
                for(int j=0;j<grid.counts[1];++j)
                    for(int k=0;k<grid.counts[2];++k)
                        aos.emplace_back(evaluator.evaluate(grid.origin+grid.steps*Eigen::Vector3d(i,j,k)));
            ao_best=std::min(ao_best,seconds(t0,Clock::now()));
        }
        std::vector<double> dense(count);
        double dense_best=1e100;
        volatile double checksum=0;
        for(int rep=0;rep<repeats;++rep) {
            const auto t0=Clock::now();
            for(std::size_t i=0;i<count;++i) dense[i]=aos[i].dot(p*aos[i]);
            dense_best=std::min(dense_best,seconds(t0,Clock::now()));
            checksum=checksum+dense[count/2];
        }
        // Experimental diagnostic only: naive AO-amplitude screening is NOT error bounded.
        std::vector<Eigen::Index> active;
        active.reserve(n);
        double screened_best=1e100, max_abs=0, sum_abs=0, sum_ref=0;
        std::size_t total_active=0;
        for(int rep=0;rep<repeats;++rep) {
            const auto t0=Clock::now();
            double local_max=0, local_sum=0, local_ref=0;
            std::size_t local_active=0;
            for(std::size_t q=0;q<count;++q) {
                active.clear();
                for(Eigen::Index i=0;i<n;++i)
                    if(std::abs(aos[q][i])>threshold) active.push_back(i);
                double rho=0;
                for(auto i:active) {
                    double row=0;
                    for(auto j:active) row+=p(i,j)*aos[q][j];
                    rho+=aos[q][i]*row;
                }
                local_active+=active.size();
                const double err=std::abs(rho-dense[q]);
                local_max=std::max(local_max,err);
                local_sum+=err;
                local_ref+=std::abs(dense[q]);
            }
            screened_best=std::min(screened_best,seconds(t0,Clock::now()));
            max_abs=local_max; sum_abs=local_sum; sum_ref=local_ref; total_active=local_active;
        }
        std::cout<<std::setprecision(12)
                 <<"basis="<<n<<" points="<<count<<" repeats="<<repeats<<" threshold="<<threshold<<"\n"
                 <<"ao_evaluation_s="<<ao_best<<"\n"
                 <<"dense_contraction_s="<<dense_best<<"\n"
                 <<"ao_plus_dense_s="<<ao_best+dense_best<<"\n"
                 <<"screened_contraction_s="<<screened_best<<"\n"
                 <<"mean_active_ao="<<double(total_active)/count<<"\n"
                 <<"screened_max_abs="<<max_abs<<"\n"
                 <<"screened_relative_L1="<<(sum_ref>0?sum_abs/sum_ref:0)<<"\n"
                 <<"checksum="<<checksum<<"\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<"\n"; return 1; }
}
