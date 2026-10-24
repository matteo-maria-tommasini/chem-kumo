#include "chem/kumo/local_neighbor_w.hpp"
#include "chem/kumo/fchk.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    if(argc<2)return 2;
    try {
        using namespace chem::kumo;
        for(int arg=1;arg<argc;++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto mol=Molecule::from_fchk(argv[arg]);
            const std::vector<Eigen::MatrixXd> densities{
                mol.density_matrix(Spin::Total),mol.density_matrix(Spin::Difference)};
            const auto full=local_density_overlap_multichannel(basis,densities);
            const auto edges=local_neighbor_pairs(mol,3.5); // bohr; test geometry only
            std::size_t checks=0;
            double worst=0;
            for(const auto& [a,b]:edges) {
                const auto block=local_neighbor_w_block(basis,densities,a,b);
                for(std::size_t ch=0;ch<densities.size();++ch)
                    for(Eigen::Index i=0;i<block.matrices[ch].rows();++i)
                        for(Eigen::Index j=0;j<block.matrices[ch].cols();++j) {
                            const double ref=full.matrices[ch](block.ao_a[static_cast<std::size_t>(i)],block.ao_b[static_cast<std::size_t>(j)]);
                            const double err=std::abs(block.matrices[ch](i,j)-ref);
                            worst=std::max(worst,err);
                            if(!std::isfinite(err)||err>1e-10*std::max(1.0,std::abs(ref)))
                                throw std::runtime_error("neighbor W mismatch");
                            ++checks;
                        }
            }
            if(edges.empty())throw std::runtime_error("no neighbor edges");
            std::cout<<"system="<<argv[arg]<<" edges="<<edges.size()<<" checks="<<checks<<" max_abs="<<worst<<" status=OK\n";
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
