// Evaluator point path across OpenMP worker counts: the blocks of grid points
// run on worker threads, and the results must be bit-identical to the
// single-thread evaluation for every block size, backend and property.
// In serial builds the same comparisons run on the calling thread.
#include "chem/kumo.hpp"
#include "omp_support.hpp"
#include "test_support.hpp"
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;

namespace {

Request full_request() {
    Request request;
    for (auto d : {density::total, density::alpha, density::beta, density::spin})
        request.add(d.value()).add(d.gradient()).add(d.hessian()).add(d.laplacian()).add(d.tau());
    request.add(orbital::alpha(1).value());
    return request;
}

// Compare every component of every property bit for bit.
void compare(const Results& expected, const Results& actual, const Request& request,
             const std::string& label) {
    kumo_test::check(expected.backend_used() == actual.backend_used(),
                     label + ": backend independent of the thread count");
    for (auto key : request.properties()) {
        const auto& e = expected.get(key);
        const auto& a = actual.get(key);
        if (!kumo_test::check(e.components.size() == a.components.size(),
                              label + ": same number of property components"))
            continue;
        std::size_t mismatches = 0;
        for (std::size_t c = 0; c < e.components.size(); ++c)
            for (std::size_t p = 0; p < e.components[c].size(); ++p)
                if (!(e.components[c][p] == a.components[c][p])) ++mismatches;
        kumo_test::check(mismatches == 0, label + ": property kind " +
                                              std::to_string(static_cast<int>(key.kind)) +
                                              " bit-identical to one thread (" +
                                              std::to_string(mismatches) + " values differ)");
    }
}

} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "file.fchk[.gz]..."))
        return 2;
    return kumo_test::run([&] {
        const Request request = full_request();
        for (int file = 1; file < argc; ++file) {
            const auto molecule = Molecule::from_fchk(argv[file]);
            const Grid grid = Grid::around(molecule.atoms(), 2.0, 0.9);
            const std::size_t n = static_cast<std::size_t>(grid.counts.prod());
            std::vector<Backend> backends{Backend::AO_P};
            if (molecule.electronic_structure_info().ao_mo_eligible())
                backends.push_back(Backend::AO_MO);
            for (auto backend : backends) {
                // Single points, uneven blocks, the default, and one block for all.
                for (std::size_t block :
                     {std::size_t(1), std::size_t(37), std::size_t(128), n + 5}) {
                    const Evaluator evaluator(molecule, {.block_size = block, .backend = backend});
                    kumo_test::set_threads(1);
                    const auto reference = evaluator.evaluate(grid, request);
                    for (int threads : {2, 4, 7}) {
                        kumo_test::set_threads(threads);
                        const std::string label = std::string(argv[file]) +
                                                  (backend == Backend::AO_MO ? " AO_MO" : " AO_P") +
                                                  " block " + std::to_string(block) + " threads " +
                                                  std::to_string(threads);
                        compare(reference, evaluator.evaluate(grid, request), request, label);
                    }
                }
            }
            // An error inside a worker block (here from the progress
            // callback, which runs in the block loop) reaches the caller.
            for (const Backend backend : {Backend::AO_P, Backend::AO_MO}) {
                if (backend == Backend::AO_MO &&
                    !molecule.electronic_structure_info().ao_mo_eligible())
                    continue;
                kumo_test::set_threads(4);
                const Evaluator evaluator(molecule, {.block_size = 7,
                                                     .backend = backend,
                                                     .progress = [](std::size_t done, std::size_t) {
                                                         if (done > 20)
                                                             throw std::domain_error("stop");
                                                     }});
                Request gradient;
                gradient.add(density::total.gradient());
                kumo_test::check_throws<std::domain_error>(
                    [&] { (void)evaluator.evaluate(grid, gradient); },
                    std::string(argv[file]) +
                        ": an error in a worker block of the point path reaches the caller");
            }
            // Block sizes outside [1, max_block_size] are rejected up front,
            // whatever the backend and the request.
            for (const std::size_t block : {std::size_t{0}, max_block_size + 1})
                kumo_test::check_throws<std::invalid_argument>(
                    [&] { (void)Evaluator(molecule, {.block_size = block}); },
                    std::string(argv[file]) + ": block_size " + std::to_string(block) +
                        " is rejected by the Evaluator");
            (void)Evaluator(molecule, {.block_size = max_block_size});
            std::cout << argv[file] << " (" << n << " points) evaluator point path OK\n";
        }
    });
}
