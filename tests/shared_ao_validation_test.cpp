#include "chem/kumo.hpp"
#include "molecule_test_access.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace chem::kumo;

namespace {
constexpr std::array<DensityComponent, 4> components{
    DensityComponent::Total, DensityComponent::Alpha, DensityComponent::Beta,
    DensityComponent::Difference};
constexpr std::array<DensityComponent, 4> spins{DensityComponent::Total, DensityComponent::Alpha,
                                                DensityComponent::Beta,
                                                DensityComponent::Difference};

// Both finite and within 5e-11 * max(1, |got|, |expected|).
void close(double got, double expected, const char* name) {
    kumo_test::check_close(got, expected, 5e-11, 5e-11, name);
}

Grid make_grid() {
    Grid grid;
    grid.origin = Eigen::Vector3d(0.13, -0.42, 0.31);
    grid.steps << 0.19, 0.02, 0.01, -0.03, 0.17, 0.02, 0.01, -0.01, 0.23;
    grid.counts = Eigen::Vector3i(5, 4, 7); // 140 points, partial block for B=128
    return grid;
}

void validate(const Molecule& molecule, std::size_t block_size) {
    const Grid grid = make_grid();
    const std::size_t n = 140;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(n));
    for (int i = 0; i < grid.counts[0]; ++i)
        for (int j = 0; j < grid.counts[1]; ++j)
            for (int k = 0; k < grid.counts[2]; ++k) {
                const auto p =
                    static_cast<std::size_t>((i * grid.counts[1] + j) * grid.counts[2] + k);
                points.col(static_cast<Eigen::Index>(p)) =
                    grid.origin + grid.steps * Eigen::Vector3d(i, j, k);
            }
    Request request;
    for (auto c : components) {
        const DensityDescriptor d{c};
        request.add(d.gradient()).add(d.hessian()).add(d.laplacian()).add(d.tau());
    }
    Evaluator evaluator(molecule,
                        {.threads = 1, .block_size = block_size, .backend = Backend::AO_P});
    const auto fused = evaluator.evaluate(grid, request);
    for (std::size_t c = 0; c < components.size(); ++c) {
        const auto reference =
            kumo_test::Kernels::density_derivatives_blocked(molecule, points, spins[c], block_size);
        const DensityDescriptor d{components[c]};
        const auto& grad = fused.get(d.gradient());
        const auto& hess = fused.get(d.hessian());
        const auto& lap = fused.get(d.laplacian());
        const auto& tau = fused.get(d.tau());
        constexpr std::array<int, 6> rows{0, 0, 0, 1, 1, 2};
        constexpr std::array<int, 6> cols{0, 1, 2, 1, 2, 2};
        for (std::size_t p = 0; p < n; ++p) {
            for (int a = 0; a < 3; ++a)
                close(grad.component(a)[p], reference[p].gradient[a],
                      "shared AO gradient matches density_derivatives_blocked()");
            for (int a = 0; a < 6; ++a)
                close(hess.component(a)[p], reference[p].hessian(rows[a], cols[a]),
                      "shared AO Hessian matches density_derivatives_blocked()");
            close(lap.component()[p], reference[p].laplacian,
                  "shared AO Laplacian matches density_derivatives_blocked()");
            close(tau.component()[p], reference[p].tau,
                  "shared AO tau matches density_derivatives_blocked()");
            close(lap.component()[p],
                  hess.component(0)[p] + hess.component(3)[p] + hess.component(5)[p],
                  "Laplacian equals the Hessian trace");
        }
    }
    // Verify spin algebra for every primitive, including the kinetic density.
    for (auto kind : {PropertyKind::Gradient, PropertyKind::Hessian, PropertyKind::Laplacian,
                      PropertyKind::Tau}) {
        const auto& total = fused.get(Property{kind, DensityComponent::Total});
        const auto& alpha = fused.get(Property{kind, DensityComponent::Alpha});
        const auto& beta = fused.get(Property{kind, DensityComponent::Beta});
        const auto& spin = fused.get(Property{kind, DensityComponent::Difference});
        for (std::size_t a = 0; a < total.components.size(); ++a)
            for (std::size_t p = 0; p < n; ++p) {
                close(total.component(a)[p], alpha.component(a)[p] + beta.component(a)[p],
                      "total = alpha + beta");
                close(spin.component(a)[p], alpha.component(a)[p] - beta.component(a)[p],
                      "spin = alpha - beta");
            }
    }
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, 2, "FILE1.fchk FILE2.fchk")) return 2;
    return kumo_test::run([&] {
        for (int i = 1; i < argc; ++i) {
            const auto molecule = Molecule::from_fchk(argv[i]);
            for (std::size_t b :
                 {std::size_t{1}, std::size_t{7}, std::size_t{128}, std::size_t{256}}) {
                const int failures_before = kumo_test::failures;
                validate(molecule, b);
                if (kumo_test::failures != failures_before)
                    std::cerr << "  (in " << argv[i] << " with block size " << b << ")\n";
            }
        }
    });
}
