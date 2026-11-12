#include "chem/kumo/evaluator.hpp"
#include "evaluation_plan.hpp"
#include "automatic_backend.hpp"
#include "ao_screening.hpp"
#include "grid_geometry.hpp"
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

Spin spin_of(SpinChannel c) {
    return c == SpinChannel::Alpha ? Spin::Alpha : Spin::Beta;
}

std::size_t channels(PropertyKind k) {
    return k == PropertyKind::Gradient ? 3 : k == PropertyKind::Hessian ? 6 : 1;
}

constexpr DensityComponent all_components[] = {DensityComponent::Total, DensityComponent::Alpha,
                                               DensityComponent::Beta, DensityComponent::Spin};

// Two evaluation paths exist:
//  - point path: shared AO derivative evaluation on explicit point blocks,
//    used for all density derivatives and, with AO-MO, also for densities;
//  - grid path: the validated Molecule::density / Molecule::orbital kernels,
//    used for AO-P densities and for all orbitals.
bool uses_point_path(const detail::ComponentPlan& component, Backend backend) {
    return component.derivatives() || (backend == Backend::AO_MO && component.density);
}

bool uses_point_path(const detail::EvaluationPlan& plan, Backend backend) {
    return std::any_of(plan.components.begin(), plan.components.end(),
                       [&](const auto& c) { return uses_point_path(c, backend); });
}

const detail::ComponentPlan& plan_for(const detail::EvaluationPlan& plan, DensityComponent c) {
    return plan.components[detail::EvaluationPlan::index(c)];
}

// Copy one property of block results, starting at grid point `start`.
// Evaluator Hessians use row-major upper-triangle order: xx, xy, xz, yy, yz, zz.
void store_block(PropertyField& output, PropertyKind kind,
                 const std::vector<DensityDerivatives>& block, std::size_t start) {
    for (std::size_t local = 0; local < block.size(); ++local) {
        const auto& d = block[local];
        const std::size_t p = start + local;
        switch (kind) {
        case PropertyKind::Density: output.components[0][p] = d.value; break;
        case PropertyKind::Gradient:
            for (int a = 0; a < 3; ++a)
                output.components[a][p] = d.gradient[a];
            break;
        case PropertyKind::Hessian:
            output.components[0][p] = d.hessian(0, 0);
            output.components[1][p] = d.hessian(0, 1);
            output.components[2][p] = d.hessian(0, 2);
            output.components[3][p] = d.hessian(1, 1);
            output.components[4][p] = d.hessian(1, 2);
            output.components[5][p] = d.hessian(2, 2);
            break;
        case PropertyKind::Laplacian: output.components[0][p] = d.laplacian; break;
        case PropertyKind::Tau: output.components[0][p] = d.tau; break;
        default: throw std::logic_error("invalid derivative property");
        }
    }
}

// Point path. Grid points are processed in blocks of at most block_size, so
// coordinate and derivative temporaries stay bounded. `evaluate_block`
// returns the derivatives of the requested spin components on one block.
template <typename BlockEvaluator>
void evaluate_point_path(const Grid& grid, std::size_t n, std::size_t block_size,
                         const detail::EvaluationPlan& plan, Backend backend,
                         const BlockEvaluator& evaluate_block,
                         std::map<Property, PropertyField>& fields) {
    std::vector<Spin> needed;
    for (auto component : all_components)
        if (uses_point_path(plan_for(plan, component), backend))
            needed.push_back(spin_of(component));

    const std::size_t block_capacity = std::min(n, block_size);
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(block_capacity));
    for (std::size_t start = 0; start < n; start += block_capacity) {
        const std::size_t count = std::min(block_capacity, n - start);
        for (std::size_t local = 0; local < count; ++local)
            points.col(static_cast<Eigen::Index>(local)) = detail::grid_point(grid, start + local);
        const auto block_points = points.leftCols(static_cast<Eigen::Index>(count)).eval();
        const auto block = evaluate_block(block_points, needed);

        for (auto component : all_components) {
            if (!uses_point_path(plan_for(plan, component), backend)) continue;
            const auto& derivatives = block.at(spin_of(component));
            for (auto& [key, output] : fields) {
                // AO-P densities come from the grid path instead.
                if (key.component != component || key.kind == PropertyKind::Orbital ||
                    (key.kind == PropertyKind::Density && backend == Backend::AO_P))
                    continue;
                store_block(output, key.kind, derivatives, start);
            }
        }
    }
}
} // namespace

Evaluator::Evaluator(const Molecule& molecule, EvaluationOptions options)
    : molecule_(&molecule), options_(options) {
    if (!options.block_size || options.threads < 0)
        throw std::invalid_argument("invalid evaluation options");
    detail::validate_screening_threshold(options.screening_threshold);
    if (options.backend != Backend::AUTO && options.backend != Backend::AO_P &&
        options.backend != Backend::AO_MO)
        throw std::invalid_argument("unknown evaluation backend");
    if (options.backend == Backend::AO_MO && !molecule.electronic_structure_info().ao_mo_eligible())
        throw std::invalid_argument(
            "AO_MO backend requested but orbital representation is not eligible");
    if (options.threads > 1)
        throw std::invalid_argument(
            "explicit multi-thread control is not implemented in the API adapter");
}

Results Evaluator::evaluate(const Grid& grid, const Request& request) const {
    // Validate the grid up front: the derivative paths below evaluate on
    // explicit points and would otherwise accept a degenerate geometry.
    const std::size_t n = detail::validated_point_count(grid);
    // Plan the requested properties before resolving the electronic backend.
    const auto plan = detail::EvaluationPlan::from(request);
    const auto& info = molecule_->electronic_structure_info();
    const Backend backend =
        detail::choose_backend(options_.backend, info.ao_mo_eligible(), molecule_->basis_size(),
                               info.occupied_alpha, info.occupied_beta, plan);
    Results result;
    result.backend_used_ = backend;
    for (auto key : request.properties()) {
        result.fields_.try_emplace(
            key, PropertyField{n, std::vector<std::vector<double>>(channels(key.kind),
                                                                   std::vector<double>(n))});
    }

    if (uses_point_path(plan, backend)) {
        // The AO-MO and AO-P block kernels are private to Molecule; Evaluator
        // is a friend, so the call is wrapped here and passed to the helper.
        const auto evaluate_block = [&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                                        const std::vector<Spin>& needed) {
            return backend == Backend::AO_MO
                       ? molecule_->density_derivatives_mo(points, needed, options_.block_size,
                                                           std::max(0, plan.maximum_ao_order),
                                                           options_.screening_threshold)
                       : molecule_->density_derivatives_shared(points, needed, options_.block_size,
                                                               plan.maximum_ao_order,
                                                               options_.screening_threshold);
        };
        evaluate_point_path(grid, n, options_.block_size, plan, backend, evaluate_block,
                            result.fields_);
    }

    // Grid path: orbitals always, densities only with AO-P.
    for (auto& [key, output] : result.fields_) {
        if (key.kind != PropertyKind::Density && key.kind != PropertyKind::Orbital) continue;
        if (key.kind == PropertyKind::Density && backend == Backend::AO_MO) continue;
        const Field field =
            key.kind == PropertyKind::Density
                ? molecule_
                      ->densities_blocked(grid, {spin_of(key.component)}, options_.block_size,
                                          options_.screening_threshold)
                      .at(spin_of(key.component))
                : molecule_->orbital(grid, key.orbital_index, spin_of(key.channel));
        std::copy(field.values().begin(), field.values().end(), output.components[0].begin());
    }
    return result;
}
} // namespace chem::kumo
