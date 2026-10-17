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
    // Build the Cartesian coordinates only when a derivative is requested.
    const bool needs_derivatives = plan.needs_derivative_points();
    Eigen::Matrix<double, 3, Eigen::Dynamic> points;
    if (needs_derivatives) {
        points.resize(3, static_cast<Eigen::Index>(n));
        for (int i = 0; i < grid.counts[0]; ++i)
            for (int j = 0; j < grid.counts[1]; ++j)
                for (int k = 0; k < grid.counts[2]; ++k) {
                    const auto idx = (static_cast<std::size_t>(i)*grid.counts[1]+j)*grid.counts[2]+k;
                    points.col(static_cast<Eigen::Index>(idx)) = grid.origin + grid.steps * Eigen::Vector3d(i,j,k);
                }
    }
    // One derivative contraction per requested spin component, independent
    // of how many of gradient/Hessian/Laplacian/tau were requested.
    for (auto component : {DensityComponent::Total, DensityComponent::Alpha,
                           DensityComponent::Beta, DensityComponent::Spin}) {
        if (!plan.components[detail::EvaluationPlan::index(component)].derivatives()) continue;
        const auto derivatives = molecule_->density_derivatives_blocked(
            points, spin_of(component), options_.block_size);
        for (auto& [key, output] : result.fields_) {
            if (key.component != component || key.kind == PropertyKind::Density ||
                key.kind == PropertyKind::Orbital) continue;
            for (std::size_t p = 0; p < n; ++p) {
                const auto& d = derivatives[p];
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
