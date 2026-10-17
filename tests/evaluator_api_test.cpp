#include "chem/kumo.hpp"
#include <cmath>
#include <stdexcept>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    using namespace chem::kumo;
    const auto molecule = Molecule::from_fchk(argv[1]);
    Grid grid;
    grid.counts = Eigen::Vector3i(2,2,2);
    grid.steps *= 0.25;
    Evaluator evaluator(molecule, {.threads=1, .block_size=4, .backend=Backend::AO_P});
    Request request;
    request.add(density::total.value()).add(density::total.gradient()).add(density::total.hessian()).add(density::total.laplacian()).add(density::total.tau());
    const auto results = evaluator.evaluate(grid, request);
    const auto legacy = molecule.density(grid);
    auto rho = results.get(density::total.value()).component();
    for (std::size_t i=0;i<rho.size();++i)
        if (std::abs(rho[i]-legacy.values()[i]) > 1e-12) throw std::runtime_error("density mismatch");
    const auto& h = results.get(density::total.hessian());
    const auto lap = results.get(density::total.laplacian()).component();
    for (std::size_t i=0;i<rho.size();++i)
        if (std::abs(lap[i]-(h.component(0)[i]+h.component(3)[i]+h.component(5)[i])) > 1e-10)
            throw std::runtime_error("Hessian trace mismatch");
}
