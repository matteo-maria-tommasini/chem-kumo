#include "chem/kumo/local_atomic_hybrids.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/local_shell_integrals.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct AOInfo { int l; std::size_t shell; };
std::vector<AOInfo> ao_metadata(const chem::kumo::BasisSet& basis) {
    using namespace chem::kumo;
    const auto components=local_shell_components(basis);
    std::vector<AOInfo> info;
    std::size_t part=0;
    for (std::size_t s=0;s<basis.shells().size();++s) {
        const auto& shell=basis.shells()[s];
        if (shell.is_sp) {
            for(int l=0;l<=1;++l) {
                for(std::size_t j=0;j<components.at(part).size();++j) info.push_back({l,s});
                ++part;
            }
        } else {
            for(std::size_t j=0;j<components.at(part).size();++j)
                info.push_back({shell.angular_momentum,s});
            ++part;
        }
    }
    if(part!=components.size() || info.size()!=basis.function_count())
        throw std::runtime_error("AO metadata ordering mismatch");
    return info;
}
const char* angular_name(int l) {
    switch(l) {case 0:return "s";case 1:return "p";case 2:return "d";case 3:return "f";case 4:return "g";default:return "higher";}
}
void print_center(const chem::kumo::BasisSet& basis,
                  const std::vector<Eigen::MatrixXd>& densities,
                  const Eigen::MatrixXd& overlap,
                  const std::vector<AOInfo>& metadata,
                  std::size_t center, std::size_t max_orbitals) {
    using namespace chem::kumo;
    const auto atomic=local_atomic_w(basis,densities,center);
    const auto hybrids=local_atomic_hybrids(atomic,overlap);
    const Eigen::Index n=hybrids.overlap.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(hybrids.overlap);
    if(eig.info()!=Eigen::Success) throw std::runtime_error("overlap eigendecomposition failed");
    const Eigen::VectorXd roots=eig.eigenvalues().cwiseMax(0.0).cwiseSqrt();
    const Eigen::MatrixXd sqrt_s=eig.eigenvectors()*roots.asDiagonal()*eig.eigenvectors().transpose();
    for(std::size_t channel=0;channel<hybrids.channels.size();++channel) {
        const auto& ch=hybrids.channels[channel];
        const Eigen::Index limit=std::min<Eigen::Index>(hybrids.rank,static_cast<Eigen::Index>(max_orbitals));
        for(Eigen::Index k=0;k<limit;++k) {
            // Symmetric Lowdin AO attribution: sum_l |(S^1/2 c)_mu|^2 = 1.
            const Eigen::VectorXd b=sqrt_s*ch.coefficients.col(k);
            double weights[6]={};
            Eigen::Index dominant=0;
            b.cwiseAbs().maxCoeff(&dominant);
            for(Eigen::Index i=0;i<n;++i) {
                const auto global=hybrids.ao_indices.at(static_cast<std::size_t>(i));
                const int l=metadata.at(static_cast<std::size_t>(global)).l;
                weights[std::clamp(l,0,5)]+=b[i]*b[i];
            }
            const auto gi=hybrids.ao_indices.at(static_cast<std::size_t>(dominant));
            const double gap_before=k?std::abs(ch.eigenvalues[k-1]-ch.eigenvalues[k]):0.;
            const double gap_after=k+1<hybrids.rank?std::abs(ch.eigenvalues[k]-ch.eigenvalues[k+1]):0.;
            const double scale=std::max(1.0,std::abs(ch.eigenvalues[k]));
            const bool near_degenerate=(k>0 && gap_before<=1e-7*scale) ||
                (k+1<hybrids.rank && gap_after<=1e-7*scale);
            std::cout<<center<<','<<channel<<','<<k<<','<<n<<','<<hybrids.rank<<','
                     <<ch.eigenvalues[k]<<','<<weights[0]<<','<<weights[1]<<','
                     <<weights[2]<<','<<weights[3]<<','<<weights[4]<<','<<weights[5]<<','
                     <<gi<<','<<metadata.at(static_cast<std::size_t>(gi)).shell<<','
                     <<angular_name(metadata.at(static_cast<std::size_t>(gi)).l)<<','
                     <<ch.coefficients(dominant,k)<<','<<b[dominant]<<','
                     <<(near_degenerate?1:0)<<','<<ch.orthogonality_max<<','<<ch.residual_max<<'\n';
        }
    }
}
} // namespace

int main(int argc,char** argv) {
    try {
        if(argc<2) {
            std::cerr<<"usage: kumo-local-atomic-hybrid-diagnostic [--all-centers] [--max-orbitals N] file.fchk[.gz]...\n";
            return 2;
        }
        bool all=false;
        std::size_t max_orbitals=8;
        std::vector<std::string> files;
        for(int i=1;i<argc;++i) {
            const std::string a=argv[i];
            if(a=="--all-centers") all=true;
            else if(a=="--max-orbitals") {
                if(++i>=argc) throw std::invalid_argument("missing --max-orbitals value");
                max_orbitals=std::stoul(argv[i]);
                if(max_orbitals==0) throw std::invalid_argument("max-orbitals must be positive");
            } else files.push_back(a);
        }
        if(files.empty()) throw std::invalid_argument("missing input files");
        std::cout<<std::setprecision(15);
        std::cout<<"center,channel,orbital,n_ao,rank,eigenvalue,weight_s,weight_p,weight_d,weight_f,weight_g,weight_higher,dominant_global_ao,dominant_shell,dominant_l,dominant_coefficient,dominant_lowdin_amplitude,near_degenerate,orthogonality_max,projected_residual_max\n";
        for(const auto& filename:files) {
            const chem::kumo::BasisSet basis{chem::kumo::FormattedCheckpoint(filename)};
            const auto molecule=chem::kumo::Molecule::from_fchk(filename);
            const std::vector<Eigen::MatrixXd> densities{
                molecule.density_matrix(chem::kumo::Spin::Total),
                molecule.density_matrix(chem::kumo::Spin::Difference)};
            const Eigen::MatrixXd overlap=chem::kumo::BasisEvaluator(basis).overlap_matrix();
            const auto metadata=ao_metadata(basis);
            if(all) {
                std::vector<std::size_t> centers;
                for(const auto& sh:basis.shells())
                    if(std::find(centers.begin(),centers.end(),sh.center)==centers.end()) centers.push_back(sh.center);
                for(const auto center:centers)
                    print_center(basis,densities,overlap,metadata,center,max_orbitals);
            } else print_center(basis,densities,overlap,metadata,basis.shells().front().center,max_orbitals);
        }
    } catch(const std::exception& e) {std::cerr<<"hybrid diagnostic: "<<e.what()<<'\n';return 1;}
}
