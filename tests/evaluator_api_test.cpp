#include "chem/kumo.hpp"
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <cmath>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, 1, "FILE.fchk")) return 2;
    using namespace chem::kumo;
    return kumo_test::run([&] {
        const auto molecule = Molecule::from_fchk(argv[1]);
        Grid grid;
        grid.counts = Eigen::Vector3i(2, 2, 2);
        grid.steps *= 0.25;
        Evaluator evaluator(molecule, {.threads = 1, .block_size = 4, .backend = Backend::AO_P});
        Request request;
        request.add(density::total.value())
            .add(density::total.gradient())
            .add(density::total.hessian())
            .add(density::total.laplacian())
            .add(density::total.tau());
        const auto results = evaluator.evaluate(grid, request);
        const auto legacy = kumo_test::Kernels::density(molecule, grid);
        auto rho = results.get(density::total.value()).component();
        // Both paths contract the same AO density matrix and differ only by
        // rounding: 1e-12 absolute is a few thousand ulps of an O(1) density.
        for (std::size_t i = 0; i < rho.size(); ++i)
            kumo_test::check(!(std::abs(rho[i] - legacy.values()[i]) > 1e-12),
                             "point " + std::to_string(i) +
                                 ": Evaluator density matches Molecule::density() within 1e-12 "
                                 "(expected " +
                                 kumo_test::num(legacy.values()[i], 17) + ", got " +
                                 kumo_test::num(rho[i], 17) + ")");
        const auto& h = results.get(density::total.hessian());
        const auto lap = results.get(density::total.laplacian()).component();
        // The trace is summed again here from the stored Hessian components,
        // so it matches the Laplacian up to rounding; second derivatives are
        // much larger than the density near the nuclei, hence the looser
        // absolute 1e-10.
        for (std::size_t i = 0; i < rho.size(); ++i) {
            const double trace = h.component(0)[i] + h.component(3)[i] + h.component(5)[i];
            kumo_test::check(!(std::abs(lap[i] - trace) > 1e-10),
                             "point " + std::to_string(i) +
                                 ": Laplacian equals the Hessian trace within 1e-10 (trace " +
                                 kumo_test::num(trace, 17) + ", Laplacian " +
                                 kumo_test::num(lap[i], 17) + ")");
        }
    });
}
