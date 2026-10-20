#include "chem/kumo/local_shell_integrals.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace chem::kumo;
using Clock = std::chrono::steady_clock;

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    try {
        const BasisSet basis{FormattedCheckpoint(argv[1])};
        const auto molecule = Molecule::from_fchk(argv[1]);
        const Eigen::MatrixXd p = molecule.density_matrix(Spin::Total);
        const auto t0 = Clock::now();
        const auto shells = local_shell_components(basis);
        const auto t1 = Clock::now();
        std::vector<Eigen::Index> offsets(shells.size()+1, 0);
        std::vector<std::size_t> primitive_terms(shells.size(), 0);
        for (std::size_t a=0; a<shells.size(); ++a) {
            offsets[a+1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
            for (const auto& component : shells[a]) primitive_terms[a] += component.size();
        }
        if (p.rows() != offsets.back() || p.cols() != offsets.back())
            throw std::runtime_error("density/basis dimension mismatch");
        using Pair = std::pair<std::size_t,std::size_t>;
        std::vector<Pair> pairs;
        for (std::size_t a=0; a<shells.size(); ++a)
            for (std::size_t b=a; b<shells.size(); ++b)
                pairs.emplace_back(a,b);
        std::size_t quartets=0, primitive_calls=0, component_integrals=0;
        double kernel_seconds=0.0, contraction_seconds=0.0, checksum=0.0;
        for (std::size_t ip=0; ip<pairs.size(); ++ip) {
            const auto [a,b] = pairs[ip];
            for (std::size_t iq=ip; iq<pairs.size(); ++iq) {
                const auto [c,d] = pairs[iq];
                const auto k0 = Clock::now();
                const auto q = local_shell_quartet({&shells[a],&shells[b],&shells[c],&shells[d]});
                const auto k1 = Clock::now();
                // Deliberately isolate contraction overhead from quartet generation.
                double block_sum = 0.0;
                for (std::size_t i=0; i<q.dimensions[0]; ++i)
                    for (std::size_t j=0; j<q.dimensions[1]; ++j)
                        for (std::size_t k=0; k<q.dimensions[2]; ++k)
                            for (std::size_t l=0; l<q.dimensions[3]; ++l)
                                block_sum += p(offsets[c]+static_cast<Eigen::Index>(k),
                                               offsets[d]+static_cast<Eigen::Index>(l)) * q.at(i,j,k,l);
                const auto k2 = Clock::now();
                checksum += block_sum;
                kernel_seconds += std::chrono::duration<double>(k1-k0).count();
                contraction_seconds += std::chrono::duration<double>(k2-k1).count();
                ++quartets;
                component_integrals += q.values.size();
                primitive_calls += primitive_terms[a]*primitive_terms[b]*primitive_terms[c]*primitive_terms[d];
            }
        }
        const double preparation_seconds=std::chrono::duration<double>(t1-t0).count();
        std::cout << std::scientific << std::setprecision(9)
                  << "system=" << argv[1] << " shells=" << shells.size()
                  << " nao=" << p.rows() << " shell_pairs=" << pairs.size()
                  << " quartets=" << quartets
                  << " component_integrals=" << component_integrals
                  << " primitive_kernel_calls=" << primitive_calls
                  << " preparation_seconds=" << preparation_seconds
                  << " quartet_seconds=" << kernel_seconds
                  << " contraction_seconds=" << contraction_seconds
                  << " checksum=" << checksum << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
