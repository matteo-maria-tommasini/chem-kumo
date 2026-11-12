// Progress callbacks of grid evaluations (densities_blocked, orbitals and
// Evaluator): calls are monotonic and end at done == total, the totals are
// the documented ones, results are bit-identical to an evaluation without a
// callback, and an exception thrown by the callback reaches the caller.
#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace chem::kumo;

// Records every call and checks that `done` increases and `total` is fixed.
struct Recorder {
    std::vector<std::size_t> done;
    std::size_t total = 0;
    bool monotonic = true;

    ProgressCallback callback() {
        return [this](std::size_t d, std::size_t t) {
            if (!done.empty() && (d <= done.back() || t != total)) monotonic = false;
            done.push_back(d);
            total = t;
        };
    }

    // Monotonic calls that end exactly at the expected total.
    void check(const std::string& label, std::size_t expected_total) const {
        kumo_test::check(!done.empty(), label + ": progress callback is called");
        kumo_test::check(monotonic, label + ": progress increases and the total is fixed");
        kumo_test::check(total == expected_total, label + ": total is " +
                                                      std::to_string(expected_total) + " (got " +
                                                      std::to_string(total) + ")");
        kumo_test::check(!done.empty() && done.back() == total,
                         label + ": the last call reports done == total");
    }
};

bool same(std::span<const double> a, std::span<const double> b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " FILE.fchk...\n";
        return 2;
    }
    for (int file = 1; file < argc; ++file) {
        const std::string name = argv[file];
        const auto molecule = Molecule::from_fchk(name);
        const auto grid = Grid::around(molecule.atoms(), 2.0, 0.6);
        const auto n = static_cast<std::size_t>(grid.counts.prod());

        // Molecule::densities_blocked: one pass over the grid.
        Recorder density_progress;
        const auto all = {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference};
        const auto densities = molecule.densities_blocked(
            grid, all, 128, default_screening_threshold, density_progress.callback());
        const auto reference = molecule.densities_blocked(grid, all);
        density_progress.check(name + " densities_blocked", n);
        for (Spin spin : all)
            kumo_test::check(same(densities.at(spin).values(), reference.at(spin).values()),
                             name + ": densities_blocked with progress is bit-identical");

        // Molecule::orbitals: one pass over the grid.
        Recorder orbital_progress;
        const std::vector<std::size_t> mos{1, 2};
        const auto orbitals =
            molecule.orbitals(grid, mos, Spin::Alpha, orbital_progress.callback());
        const auto orbital_reference = molecule.orbitals(grid, mos, Spin::Alpha);
        orbital_progress.check(name + " orbitals", n);
        for (std::size_t i = 0; i < mos.size(); ++i)
            kumo_test::check(same(orbitals[i].values(), orbital_reference[i].values()),
                             name + ": orbitals with progress are bit-identical");

        // Evaluator: AO_P runs the derivative pass plus one grid pass per
        // density field (here 2); AO_MO evaluates everything in one pass.
        Request request;
        request.add(density::total.value())
            .add(density::spin.value())
            .add(density::total.laplacian())
            .add(density::alpha.tau());
        for (Backend backend : {Backend::AO_P, Backend::AO_MO}) {
            if (backend == Backend::AO_MO && !molecule.electronic_structure_info().ao_mo_eligible())
                continue;
            const std::string label =
                name + (backend == Backend::AO_P ? " Evaluator AO_P" : " Evaluator AO_MO");
            Recorder evaluator_progress;
            EvaluationOptions options;
            options.backend = backend;
            options.progress = evaluator_progress.callback();
            const auto results = Evaluator(molecule, options).evaluate(grid, request);
            options.progress = {};
            const auto plain = Evaluator(molecule, options).evaluate(grid, request);
            evaluator_progress.check(label, (backend == Backend::AO_P ? 3 : 1) * n);
            bool identical = true;
            for (const auto& property : request.properties())
                for (std::size_t c = 0; c < results.get(property).components.size(); ++c)
                    identical = identical && same(results.get(property).component(c),
                                                  plain.get(property).component(c));
            kumo_test::check(identical, label + ": results with progress are bit-identical");
        }

        // An exception thrown by the callback aborts the evaluation.
        const ProgressCallback stop = [](std::size_t done, std::size_t) {
            if (done > 0) throw std::runtime_error("stop");
        };
        kumo_test::check_throws<std::runtime_error>(
            [&] { (void)molecule.densities_blocked(grid, {Spin::Total}, 128, 0.0, stop); },
            name + ": a throwing progress callback aborts densities_blocked");
        kumo_test::check_throws<std::runtime_error>(
            [&] { (void)molecule.orbitals(grid, {1}, Spin::Alpha, stop); },
            name + ": a throwing progress callback aborts orbitals");
        kumo_test::check_throws<std::runtime_error>(
            [&] {
                EvaluationOptions options;
                options.progress = stop;
                (void)Evaluator(molecule, options).evaluate(grid, request);
            },
            name + ": a throwing progress callback aborts Evaluator::evaluate");
    }
    return kumo_test::exit_code();
}
