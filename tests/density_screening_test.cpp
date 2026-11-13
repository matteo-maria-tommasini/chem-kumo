// AO screening of grid densities and density derivatives (AO_P and AO_MO
// contractions): with the
// default threshold every result matches the exact contraction (threshold 0)
// to rounding level, a large threshold really changes the result (so
// screening is exercised), and invalid thresholds are rejected.
#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
double max_difference(std::span<const double> a, std::span<const double> b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        result = std::max(result, std::abs(a[i] - b[i]));
    return result;
}

double max_magnitude(std::span<const double> a) {
    double result = 0.0;
    for (double value : a)
        result = std::max(result, std::abs(value));
    return result;
}

// All numbers of one DensityDerivatives, flattened for comparisons.
std::vector<double> flatten(const chem::kumo::DensityDerivatives& d) {
    std::vector<double> out{d.value, d.laplacian, d.tau};
    for (int i = 0; i < 3; ++i)
        out.push_back(d.gradient[i]);
    for (int i = 0; i < 9; ++i)
        out.push_back(d.hessian(i));
    return out;
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        using namespace chem::kumo;
        const auto all = {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference};
        const char* names[] = {"total", "alpha", "beta", "difference"};
        for (int n = 1; n < argc; ++n) {
            const std::string file = argv[n];
            const auto molecule = Molecule::from_fchk(file);
            const auto grid = Grid::around(molecule.atoms(), 3.0, 0.8);
            const auto exact = molecule.densities_blocked(grid, all, 128, 0.0);
            const auto screened = molecule.densities_blocked(grid, all);
            const auto coarse = molecule.densities_blocked(grid, all, 128, 1e-3);
            double coarse_change = 0.0;
            for (Spin spin : all) {
                const auto label = file + " " + names[static_cast<int>(spin)];
                const double error =
                    max_difference(screened.at(spin).values(), exact.at(spin).values());
                kumo_test::check(error <= 1e-10, label +
                                                     ": default screening matches the exact "
                                                     "density within 1e-10 (max error " +
                                                     kumo_test::num(error) + ")");
                coarse_change = std::max(coarse_change, max_difference(coarse.at(spin).values(),
                                                                       exact.at(spin).values()));
            }
            kumo_test::check(coarse_change > 1e-8,
                             file +
                                 ": threshold 1e-3 changes the density, so AOs are screened "
                                 "(max change " +
                                 kumo_test::num(coarse_change) + ")");

            // Evaluator AO_P densities use the same screening and option.
            Request request;
            request.add(density::total.value()).add(density::spin.value());
            const auto evaluate = [&](double threshold) {
                EvaluationOptions options;
                options.backend = Backend::AO_P;
                options.screening_threshold = threshold;
                return Evaluator(molecule, options).evaluate(grid, request);
            };
            const auto evaluated = evaluate(default_screening_threshold);
            const auto evaluated_exact = evaluate(0.0);
            for (const auto& property : {density::total.value(), density::spin.value()}) {
                const double error = max_difference(evaluated.get(property).component(),
                                                    evaluated_exact.get(property).component());
                kumo_test::check(error <= 1e-10, file +
                                                     ": Evaluator AO_P screening matches the exact "
                                                     "density within 1e-10 (max error " +
                                                     kumo_test::num(error) + ")");
            }
            const double exact_match = max_difference(
                evaluated_exact.get(density::total.value()).component(),
                molecule.densities_blocked(grid, {Spin::Total}, 128, 0.0).at(Spin::Total).values());
            kumo_test::check(exact_match == 0.0,
                             file +
                                 ": Evaluator with threshold 0 equals densities_blocked with "
                                 "threshold 0 (max difference " +
                                 kumo_test::num(exact_match) + ")");

            // Evaluator density derivatives, with each backend available for this
            // molecule (AO_MO also contracts densities on this point path). The
            // scale of derivatives grows near the nuclei, so the error is
            // measured relative to the largest value of each component.
            Request derivatives;
            const Property derivative_properties[] = {
                density::total.value(),   density::spin.value(),      density::total.gradient(),
                density::total.hessian(), density::total.laplacian(), density::total.tau(),
                density::spin.gradient(), density::alpha.hessian()};
            for (const auto& property : derivative_properties)
                derivatives.add(property);
            std::vector<Backend> backends{Backend::AO_P};
            if (molecule.electronic_structure_info().ao_mo_eligible())
                backends.push_back(Backend::AO_MO);
            for (Backend backend : backends) {
                const std::string backend_name = backend == Backend::AO_P ? "AO_P" : "AO_MO";
                const auto evaluate_derivatives = [&](double threshold) {
                    EvaluationOptions options;
                    options.backend = backend;
                    options.screening_threshold = threshold;
                    return Evaluator(molecule, options).evaluate(grid, derivatives);
                };
                const auto derivatives_screened = evaluate_derivatives(default_screening_threshold);
                const auto derivatives_exact = evaluate_derivatives(0.0);
                const auto derivatives_coarse = evaluate_derivatives(1e-3);
                double derivative_error = 0.0, derivative_change = 0.0;
                for (const auto& property : derivative_properties) {
                    const auto& exact_field = derivatives_exact.get(property);
                    for (std::size_t c = 0; c < exact_field.components.size(); ++c) {
                        const auto exact_values = exact_field.component(c);
                        const double scale = std::max(1.0, max_magnitude(exact_values));
                        derivative_error =
                            std::max(derivative_error,
                                     max_difference(derivatives_screened.get(property).component(c),
                                                    exact_values) /
                                         scale);
                        derivative_change =
                            std::max(derivative_change,
                                     max_difference(derivatives_coarse.get(property).component(c),
                                                    exact_values) /
                                         scale);
                    }
                }
                kumo_test::check(derivative_error <= 1e-12,
                                 file + ": Evaluator " + backend_name +
                                     " screening matches the exact contraction within 1e-12 "
                                     "relative (max " +
                                     kumo_test::num(derivative_error) + ")");
                kumo_test::check(derivative_change > 1e-10,
                                 file + ": threshold 1e-3 changes the " + backend_name +
                                     " results, so AOs are screened (max relative change " +
                                     kumo_test::num(derivative_change) + ")");
            }

            // Molecule::density_derivatives_blocked uses the same screening; test
            // it on 300 points along the diagonal of the grid box.
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 300);
            const Eigen::Vector3d diagonal =
                grid.steps * (grid.counts.array() - 1).cast<double>().matrix();
            for (Eigen::Index j = 0; j < points.cols(); ++j)
                points.col(j) = grid.origin + diagonal * (static_cast<double>(j) / 299.0);
            const auto blocked_exact =
                molecule.density_derivatives_blocked(points, Spin::Total, 128, 0.0);
            const auto blocked_screened = molecule.density_derivatives_blocked(points, Spin::Total);
            double blocked_error = 0.0, blocked_scale = 1.0;
            for (std::size_t j = 0; j < blocked_exact.size(); ++j) {
                const auto a = flatten(blocked_exact[j]), b = flatten(blocked_screened[j]);
                blocked_error = std::max(blocked_error, max_difference(a, b));
                blocked_scale = std::max(blocked_scale, max_magnitude(a));
            }
            kumo_test::check(blocked_error <= 1e-12 * blocked_scale,
                             file +
                                 ": density_derivatives_blocked screening matches the exact "
                                 "contraction within 1e-12 relative (max " +
                                 kumo_test::num(blocked_error / blocked_scale) + ")");

            for (double bad : {-1e-12, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
                kumo_test::check_throws<std::invalid_argument>(
                    [&] { (void)molecule.densities_blocked(grid, {Spin::Total}, 128, bad); },
                    "densities_blocked rejects screening threshold " + kumo_test::num(bad));
                kumo_test::check_throws<std::invalid_argument>(
                    [&] {
                        (void)molecule.density_derivatives_blocked(points, Spin::Total, 128, bad);
                    },
                    "density_derivatives_blocked rejects screening threshold " +
                        kumo_test::num(bad));
                kumo_test::check_throws<std::invalid_argument>(
                    [&] {
                        EvaluationOptions options;
                        options.screening_threshold = bad;
                        Evaluator evaluator(molecule, options);
                    },
                    "Evaluator rejects screening threshold " + kumo_test::num(bad));
            }
        }
    });
}
