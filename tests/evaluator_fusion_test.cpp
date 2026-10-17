#include "chem/kumo.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;

namespace {
void require_close(double actual, double expected, const std::string& label) {
    const double tolerance = 2e-11 * std::max({1.0, std::abs(actual), std::abs(expected)});
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(label + ": mismatch: " + std::to_string(actual) + " versus " + std::to_string(expected));
}
Spin legacy_spin(DensityComponent c) {
    switch (c) {
        case DensityComponent::Total: return Spin::Total;
        case DensityComponent::Alpha: return Spin::Alpha;
        case DensityComponent::Beta: return Spin::Beta;
        case DensityComponent::Spin: return Spin::Difference;
    }
    throw std::logic_error("invalid spin component");
}
std::vector<Property> keys(DensityComponent c) {
    const DensityDescriptor d{c};
    return {d.value(), d.gradient(), d.hessian(), d.laplacian(), d.tau()};
}
void check(const Molecule& molecule, std::size_t block_size, bool reversed) {
    Grid grid;
    grid.origin = Eigen::Vector3d(0.17, -0.28, 0.36);
    grid.steps << 0.33, 0.04, 0.00,
                 -0.02, 0.29, 0.03,
                  0.01, 0.00, 0.27;
    grid.counts = Eigen::Vector3i(3, 2, 4);
    const std::size_t n = 24;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(n));
    for (int i = 0; i < grid.counts[0]; ++i)
        for (int j = 0; j < grid.counts[1]; ++j)
            for (int k = 0; k < grid.counts[2]; ++k) {
                const std::size_t p = (static_cast<std::size_t>(i)*grid.counts[1]+j)*grid.counts[2]+k;
                points.col(static_cast<Eigen::Index>(p)) = grid.origin + grid.steps * Eigen::Vector3d(i,j,k);
            }
    Evaluator evaluator(molecule, {.threads=1, .block_size=block_size, .backend=Backend::AO_P});
    Request combined;
    std::vector<Property> all;
    for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                           DensityComponent::Beta, DensityComponent::Spin}) {
        const auto props = keys(component);
        all.insert(all.end(), props.begin(), props.end());
    }
    if (reversed) std::reverse(all.begin(), all.end());
    for (const auto key : all) combined.add(key);
    const auto fused = evaluator.evaluate(grid, combined);
    for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                           DensityComponent::Beta, DensityComponent::Spin}) {
        const auto reference = molecule.density_derivatives_blocked(points, legacy_spin(component), block_size);
        const auto scalar = molecule.density(grid, legacy_spin(component));
        for (const auto key : keys(component)) {
            Request individual;
            individual.add(key);
            const auto single = evaluator.evaluate(grid, individual);
            const auto& field = fused.get(key);
            const auto& separate = single.get(key);
            if (field.points != n || field.components.size() != separate.components.size())
                throw std::runtime_error("incorrect field shape");
            for (std::size_t axis = 0; axis < field.components.size(); ++axis)
                for (std::size_t p = 0; p < n; ++p)
                    require_close(field.components[axis][p], separate.components[axis][p], "fused vs individual");
            for (std::size_t p = 0; p < n; ++p) {
                const auto& d = reference[p];
                switch (key.kind) {
                    case PropertyKind::Density:
                        require_close(field.component()[p], scalar.values()[p], "scalar density");
                        break;
                    case PropertyKind::Gradient:
                        for (int a = 0; a < 3; ++a)
                            require_close(field.component(a)[p], d.gradient[a], "gradient");
                        break;
                    case PropertyKind::Hessian:
                        for (int a = 0; a < 6; ++a) {
                            constexpr std::array<int,6> row{0,0,0,1,1,2};
                            constexpr std::array<int,6> col{0,1,2,1,2,2};
                            require_close(field.component(a)[p], d.hessian(row[a], col[a]), "Hessian");
                        }
                        break;
                    case PropertyKind::Laplacian:
                        require_close(field.component()[p], d.laplacian, "Laplacian");
                        break;
                    case PropertyKind::Tau:
                        require_close(field.component()[p], d.tau, "tau");
                        break;
                    default: throw std::logic_error("unexpected property");
                }
            }
        }
        const auto& h = fused.get(DensityDescriptor{component}.hessian());
        const auto& lap = fused.get(DensityDescriptor{component}.laplacian());
        for (std::size_t p = 0; p < n; ++p)
            require_close(lap.component()[p], h.component(0)[p] + h.component(3)[p] + h.component(5)[p], "Hessian trace");
    }
    // Duplicate descriptors must not duplicate result fields or change values.
    Request duplicate;
    duplicate.add(density::total.gradient()).add(density::total.gradient());
    const auto repeated = evaluator.evaluate(grid, duplicate);
    for (std::size_t p = 0; p < n; ++p)
        for (int a = 0; a < 3; ++a)
            require_close(repeated.get(density::total.gradient()).component(a)[p],
                          fused.get(density::total.gradient()).component(a)[p], "duplicate descriptor");
}
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    try {
        for (int file = 1; file < argc; ++file) {
            const auto molecule = Molecule::from_fchk(argv[file]);
            for (std::size_t block_size : {std::size_t{1}, std::size_t{5}, std::size_t{128}})
                for (bool reversed : {false, true}) check(molecule, block_size, reversed);
        }
    } catch (const std::exception& e) {
        std::cerr << "evaluator fusion regression: " << e.what() << '\n';
        return 1;
    }
}
