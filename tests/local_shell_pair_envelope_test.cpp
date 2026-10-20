#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/fchk.hpp"
#include "chem/kumo/basis.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    using namespace chem::kumo;
    if (argc < 2) return 2;
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            const std::size_t n = std::min<std::size_t>(shells.size(), 6);
            std::vector<LocalShellPairEnvelope> pairs;
            pairs.reserve(n*n);
            for (std::size_t a=0; a<n; ++a)
                for (std::size_t b=0; b<n; ++b)
                    pairs.push_back(prepare_local_shell_pair_envelope(shells[a],shells[b]));
            std::size_t quartets=0, elements=0;
            double violation=0, worst_ratio=0;
            for (std::size_t a=0; a<n; ++a)
                for (std::size_t b=0; b<n; ++b)
                    for (std::size_t c=0; c<n; ++c)
                        for (std::size_t d=0; d<n; ++d) {
                            const auto q=local_shell_quartet({&shells[a],&shells[b],&shells[c],&shells[d]});
                            const auto& ab=pairs[a*n+b];
                            const auto& cd=pairs[c*n+d];
                            ++quartets;
                            for (std::size_t i=0;i<ab.na;++i)
                                for (std::size_t j=0;j<ab.nb;++j)
                                    for (std::size_t k=0;k<cd.na;++k)
                                        for (std::size_t l=0;l<cd.nb;++l) {
                                            const double exact=std::abs(q.at(i,j,k,l));
                                            const double bound=local_shell_pair_envelope_bound(ab,i,j,cd,k,l);
                                            if (!std::isfinite(bound)) throw std::runtime_error("nonfinite bound");
                                            violation=std::max(violation,exact-bound);
                                            if (bound+2e-11*std::max(1.0,exact)<exact)
                                                throw std::runtime_error("envelope violation");
                                            if (exact>1e-12) worst_ratio=std::max(worst_ratio,bound/exact);
                                            ++elements;
                                        }
                        }
            std::cout<<std::setprecision(12)<<"system="<<argv[arg]
                <<" shells_tested="<<n<<" quartets="<<quartets
                <<" elements="<<elements<<" max_positive_violation="<<violation
                <<" max_bound_to_integral="<<worst_ratio<<'\n';
        }
    } catch (const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
