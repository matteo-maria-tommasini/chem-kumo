// Evaluator point path across OpenMP worker counts: the blocks of grid points
// run on worker threads, and the results must be bit-identical to the
// single-thread evaluation for every block size, backend and property.
// In serial builds the same comparisons run on the calling thread.
#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace chem::kumo;

namespace {

void set_threads(int threads) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
}

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
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " file.fchk[.gz]...\n";
        return 2;
    }
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
            for (std::size_t block : {std::size_t(1), std::size_t(37), std::size_t(128), n + 5}) {
                const Evaluator evaluator(molecule, {.block_size = block, .backend = backend});
                set_threads(1);
                const auto reference = evaluator.evaluate(grid, request);
                for (int threads : {2, 4, 7}) {
                    set_threads(threads);
                    const std::string label =
                        std::string(argv[file]) + (backend == Backend::AO_MO ? " AO_MO" : " AO_P") +
                        " block " + std::to_string(block) + " threads " + std::to_string(threads);
                    compare(reference, evaluator.evaluate(grid, request), request, label);
                }
            }
        }
        // A kernel error inside a worker block reaches the caller.
        if (molecule.electronic_structure_info().ao_mo_eligible()) {
            set_threads(4);
            const Evaluator evaluator(molecule, {.block_size = 70000, .backend = Backend::AO_MO});
            Request gradient;
            gradient.add(density::total.gradient());
            kumo_test::check_throws<std::invalid_argument>(
                [&] { (void)evaluator.evaluate(grid, gradient); },
                std::string(argv[file]) + ": invalid AO-MO block size throws from the point path");
        }
        std::cout << argv[file] << " (" << n << " points) evaluator point path OK\n";
    }
    return kumo_test::exit_code();
}
