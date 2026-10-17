#include "chem/kumo/evaluator.hpp"
#include "evaluation_plan.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace chem::kumo {
namespace {
Spin spin_of(DensityComponent c) {
    switch (c) {
        case DensityComponent::Total: return Spin::Total;
        case DensityComponent::Alpha: return Spin::Alpha;
        case DensityComponent::Beta: return Spin::Beta;
        case DensityComponent::Spin: return Spin::Difference;
    }
    throw std::invalid_argument("unknown density component");
}
std::size_t channels(PropertyKind k) {
    return k == PropertyKind::Gradient ? 3 : k == PropertyKind::Hessian ? 6 : 1;
}
}
Evaluator::Evaluator(const Molecule& molecule, EvaluationOptions options)
    : molecule_(&molecule), options_(options) {
    if (!options.block_size || options.threads < 0) throw std::invalid_argument("invalid evaluation options");
if (options.threads > 1) throw std::invalid_argument("explicit multi-thread control is not implemented in the API adapter");
}
Results Evaluator::evaluate(const Grid& grid, const Request& request) const {
    // Stage 02A: correctness-first adapter to the validated 1.3 kernels.
    // It does not yet implement fused multi-property contraction or AO-MO selection.
    std::size_t n = 1;
    for (int axis = 0; axis < 3; ++axis) {
        if (grid.counts[axis] <= 0 || n > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(grid.counts[axis]))
            throw std::invalid_argument("invalid grid size");
        n *= static_cast<std::size_t>(grid.counts[axis]);
    }
    const auto plan = detail::EvaluationPlan::from(request);
    Results result;
    for (auto key : request.properties()) {
        result.fields_.try_emplace(key, PropertyField{n, std::vector<std::vector<double>>(channels(key.kind), std::vector<double>(n))});
    }
    // Stage 02C.2a: bound coordinate and derivative temporaries by block_size.
    // The validated AO-P kernel remains the numerical reference. This is not
    // yet a fused AO evaluation across electronic components.
    if (plan.needs_derivative_points()) {
        const std::size_t block_capacity = std::min(n, options_.block_size);
        Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(block_capacity));
        const std::size_t ny = static_cast<std::size_t>(grid.counts[1]);
        const std::size_t nz = static_cast<std::size_t>(grid.counts[2]);
        for (std::size_t start = 0; start < n; start += block_capacity) {
            const std::size_t count = std::min(block_capacity, n - start);
            for (std::size_t local = 0; local < count; ++local) {
                const std::size_t linear = start + local;
                const std::size_t i = linear / (ny * nz);
                const std::size_t j = (linear / nz) % ny;
                const std::size_t k = linear % nz;
                points.col(static_cast<Eigen::Index>(local)) = grid.origin + grid.steps *
                    Eigen::Vector3d(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
            }
            const auto block_points = points.leftCols(static_cast<Eigen::Index>(count)).eval();
            std::vector<Spin> needed;
            for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                                   DensityComponent::Beta, DensityComponent::Spin}) {
                if (plan.components[detail::EvaluationPlan::index(component)].derivatives())
                    needed.push_back(spin_of(component));
            }
            const auto shared = molecule_->density_derivatives_shared(block_points, needed, options_.block_size);
            for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                                   DensityComponent::Beta, DensityComponent::Spin}) {
                if (!plan.components[detail::EvaluationPlan::index(component)].derivatives()) continue;
                const auto& derivatives = shared.at(spin_of(component));
                for (auto& [key, output] : result.fields_) {
                    if (key.component != component || key.kind == PropertyKind::Density ||
                        key.kind == PropertyKind::Orbital) continue;
                    for (std::size_t local = 0; local < count; ++local) {
                        const auto& d = derivatives[local];
                        const std::size_t p = start + local;
                        switch (key.kind) {
                            case PropertyKind::Gradient:
                                for (int a=0;a<3;++a) output.components[a][p]=d.gradient[a];
                                break;
                            case PropertyKind::Hessian:
                                output.components[0][p]=d.hessian(0,0);
                                output.components[1][p]=d.hessian(0,1);
                                output.components[2][p]=d.hessian(0,2);
                                output.components[3][p]=d.hessian(1,1);
                                output.components[4][p]=d.hessian(1,2);
                                output.components[5][p]=d.hessian(2,2);
                                break;
                            case PropertyKind::Laplacian: output.components[0][p]=d.laplacian; break;
                            case PropertyKind::Tau: output.components[0][p]=d.tau; break;
                            default: throw std::logic_error("invalid derivative property");
                        }
                    }
                }
            }
        }
    }
    // Scalar density and orbital paths remain delegated to the 1.3 kernels.
    for (auto& [key, output] : result.fields_) {
        if (key.kind != PropertyKind::Density && key.kind != PropertyKind::Orbital) continue;
        Field field = key.kind == PropertyKind::Density
            ? molecule_->density(grid, spin_of(key.component))
            : molecule_->orbital(grid, key.orbital_index,
                key.channel == SpinChannel::Alpha ? Spin::Alpha : Spin::Beta);
        std::copy(field.values().begin(), field.values().end(), output.components[0].begin());
    }
    return result;
}
}
