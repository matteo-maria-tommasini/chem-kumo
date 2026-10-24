#include "chem/kumo/local_atomic_w.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        using namespace chem::kumo;
        for(int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto molecule=Molecule::from_fchk(argv[arg]);
            const std::vector<Eigen::MatrixXd> densities{
                molecule.density_matrix(Spin::Total),
                molecule.density_matrix(Spin::Difference)};
            const auto full=local_density_overlap_multichannel(basis,densities);
            std::vector<std::size_t> centers;
            for(const auto& s:basis.shells())
                if(std::find(centers.begin(),centers.end(),s.center)==centers.end())
                    centers.push_back(s.center);
            std::size_t checks=0;
            double max_error=0.0;
            for(const auto center:centers) {
                const auto local=local_atomic_w(basis,densities,center);
                if(local.matrices.size()!=densities.size() || local.ao_indices.empty())
                    throw std::runtime_error("atomic W: missing block");
                for(std::size_t ch=0;ch<densities.size();++ch)
                    for(Eigen::Index i=0;i<local.matrices[ch].rows();++i)
                        for(Eigen::Index j=0;j<local.matrices[ch].cols();++j) {
                            const double ref=full.matrices[ch](local.ao_indices[static_cast<std::size_t>(i)],
                                                                local.ao_indices[static_cast<std::size_t>(j)]);
                            const double got=local.matrices[ch](i,j);
                            const double error=std::abs(ref-got);
                            max_error=std::max(max_error,error);
                            if(!std::isfinite(got) || error>1e-10*std::max(1.0,std::abs(ref)))
                                throw std::runtime_error("atomic W: reference mismatch");
                            ++checks;
                        }
            }
            std::cout<<"system="<<argv[arg]<<" atoms="<<centers.size()
                     <<" checks="<<checks<<" max_abs="<<max_error<<" status=OK\n";
        }
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
