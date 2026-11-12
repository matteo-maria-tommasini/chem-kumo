#include "chem/kumo/evaluator.hpp"
#include "evaluation_plan.hpp"
#include "automatic_backend.hpp"
#include "electrostatic_potential.hpp"
#include "ao_screening.hpp"
#include "grid_geometry.hpp"
#include "parallel_blocks.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>
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

// Electron localization function (Becke and Edgecombe 1990) with
// tau = (1/2) sum_i |grad psi_i|^2. Each spin channel s with rho_s > 0
// contributes the Pauli kinetic energy density
//   D_s = tau_s - |grad rho_s|^2 / (8 rho_s)
// and the Thomas-Fermi reference D0_s = (3/10) (6 pi^2)^(2/3) rho_s^(5/3);
// ELF = 1 / (1 + (D / D0)^2) with D and D0 summed over the channels given.
// For a closed shell this is the familiar spin-unpolarized form
// D = tau - |grad rho|^2 / (8 rho), D0 = (3/10) (3 pi^2)^(2/3) rho^(5/3).
// Channels without density contribute nothing; with no density at all the
// ELF is 0. The result is finite and in [0, 1] for any finite input.
double elf(std::initializer_list<const DensityDerivatives*> channels) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    static const double fermi = 0.3 * std::pow(6.0 * pi * pi, 2.0 / 3.0);
    double pauli = 0.0;
    double uniform = 0.0;
    for (const DensityDerivatives* d : channels) {
        if (!(d->value > 0.0)) continue;
        pauli += d->tau - d->gradient.squaredNorm() / (8.0 * d->value);
        uniform += fermi * std::pow(d->value, 5.0 / 3.0);
    }
    if (!(uniform > 0.0)) return 0.0;
    const double chi = pauli / uniform;
    return 1.0 / (1.0 + chi * chi);
}

// Reduced density gradient s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)),
// capped at reduced_gradient_cap (also where rho <= 0 or rho^(4/3)
// underflows, where s is undefined or infinite).
double reduced_gradient(const DensityDerivatives& d) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    static const double scale = 2.0 * std::cbrt(3.0 * pi * pi);
    if (!(d.value > 0.0)) return reduced_gradient_cap;
    const double denominator = scale * std::pow(d.value, 4.0 / 3.0);
    const double norm = d.gradient.norm();
    if (!(denominator > 0.0) || !(norm < reduced_gradient_cap * denominator))
        return reduced_gradient_cap;
    return norm / denominator;
}

// sign(lambda2) * rho, with lambda1 <= lambda2 <= lambda3 the eigenvalues of
// the density Hessian; lambda2 == 0 counts as positive.
double sign_lambda2_density(const DensityDerivatives& d) {
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(d.hessian, Eigen::EigenvaluesOnly);
    return solver.eigenvalues()[1] < 0.0 ? -d.value : d.value;
}

using BlockDerivatives = std::map<Spin, std::vector<DensityDerivatives>>;

// Copy (or derive) one property of block results, starting at grid point
// `start`. Evaluator Hessians use row-major upper-triangle order: xx, xy, xz,
// yy, yz, zz.
void store_block(PropertyField& output, const Property& key, const BlockDerivatives& block,
                 std::size_t count, std::size_t start) {
    if (key.kind == PropertyKind::Elf && key.component == DensityComponent::Total) {
        const auto& alpha = block.at(Spin::Alpha);
        const auto& beta = block.at(Spin::Beta);
        for (std::size_t local = 0; local < count; ++local)
            output.components[0][start + local] = elf({&alpha[local], &beta[local]});
        return;
    }
    const auto& derivatives = block.at(spin_of(key.component));
    for (std::size_t local = 0; local < count; ++local) {
        const auto& d = derivatives[local];
        const std::size_t p = start + local;
        switch (key.kind) {
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
        case PropertyKind::Elf: output.components[0][p] = elf({&d}); break;
        case PropertyKind::ReducedGradient: output.components[0][p] = reduced_gradient(d); break;
        case PropertyKind::SignLambda2Density:
            output.components[0][p] = sign_lambda2_density(d);
            break;
        default: throw std::logic_error("invalid derivative property");
        }
    }
}

// Point path. Grid points are processed in blocks of at most block_size, so
// coordinate and derivative temporaries stay bounded. `evaluate_block`
// returns the derivatives of the requested spin components on one block.
//
// Each kernel call receives exactly one block, so the block boundaries and
// the arithmetic per block do not depend on the number of threads. With
// OpenMP the blocks run on worker threads (the kernels then stay serial,
// since their own parallel regions open only for more than one block); each
// block writes a disjoint range of grid points, so the output is
// bit-identical to the serial loop.
template <typename BlockEvaluator>
void evaluate_point_path(const Grid& grid, std::size_t n, std::size_t block_size,
                         const detail::EvaluationPlan& plan, Backend backend,
                         const BlockEvaluator& evaluate_block,
                         std::map<Property, PropertyField>& fields,
                         const ProgressCallback& progress) {
    std::vector<Spin> needed;
    for (auto component : all_components)
        if (uses_point_path(plan_for(plan, component), backend))
            needed.push_back(spin_of(component));

    const std::size_t nblocks =
        detail::parallel_block_count(n, block_size, "too many point blocks for evaluation");
    detail::FirstException worker_errors;
    detail::Progress report(progress, n);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
#endif
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(nblocks); ++index) {
        try {
            const std::size_t start = static_cast<std::size_t>(index) * block_size;
            const std::size_t count = std::min(block_size, n - start);
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(count));
            for (std::size_t local = 0; local < count; ++local)
                points.col(static_cast<Eigen::Index>(local)) =
                    detail::grid_point(grid, start + local);
            const auto block = evaluate_block(points, needed);

            for (auto& [key, output] : fields) {
                // Orbitals and AO-P densities come from the grid path, and
                // the electrostatic potential from its own pass.
                if (key.kind == PropertyKind::Orbital ||
                    key.kind == PropertyKind::ElectrostaticPotential ||
                    (key.kind == PropertyKind::Density && backend == Backend::AO_P))
                    continue;
                store_block(output, key, block, count, start);
            }
            // Serialized and monotonic; a throwing callback is captured
            // like any other block error.
            report.advance(count);
        } catch (...) {
            worker_errors.capture();
        }
    }
    worker_errors.rethrow_if_any();
}

// Electrostatic potential pass, in blocks of grid points scheduled like the
// point path; each point is computed independently with a fixed summation
// order, so the values do not depend on the number of threads.
// Electrostatic potential pass. The grid is split into clusters of up to
// cluster_side^3 points (index boxes), each evaluated around its own centre
// with one far-field plan whose radius covers the largest cluster. Cluster
// boundaries depend only on the grid, and each value only on its point, its
// cluster and the plan, so the result does not depend on the number of
// threads or on the block size.
constexpr int cluster_side = 8;

void evaluate_electrostatic_potential(const Grid& grid, std::size_t n, double tolerance,
                                      const detail::ElectrostaticPotential& potential,
                                      std::vector<double>& output,
                                      const ProgressCallback& progress) {
    Eigen::Vector3i tiles, extent;
    for (int a = 0; a < 3; ++a) {
        extent[a] = std::min(cluster_side, grid.counts[a]);
        tiles[a] = (grid.counts[a] + cluster_side - 1) / cluster_side;
    }
    // Largest distance from a cluster centre to its points: a corner of the
    // largest index box.
    double radius = 0.0;
    for (int corner = 0; corner < 8; ++corner) {
        Eigen::Vector3d offset;
        for (int a = 0; a < 3; ++a)
            offset[a] = ((corner >> a) & 1 ? 0.5 : -0.5) * (extent[a] - 1);
        radius = std::max(radius, (grid.steps * offset).norm());
    }
    // Slightly enlarged to cover rounding in the point coordinates.
    const auto plan = potential.far_field_plan(radius * (1.0 + 1e-10) + 1e-12, tolerance);

    const auto ntiles = static_cast<std::int64_t>(tiles[0]) * tiles[1] * tiles[2];
    detail::FirstException worker_errors;
    detail::Progress report(progress, n);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (ntiles > 1)
#endif
    for (std::int64_t index = 0; index < ntiles; ++index) {
        try {
            const Eigen::Vector3i tile(static_cast<int>(index / (tiles[1] * tiles[2])),
                                       static_cast<int>((index / tiles[2]) % tiles[1]),
                                       static_cast<int>(index % tiles[2]));
            Eigen::Vector3i begin, end;
            for (int a = 0; a < 3; ++a) {
                begin[a] = tile[a] * cluster_side;
                end[a] = std::min(begin[a] + cluster_side, grid.counts[a]);
            }
            const Eigen::Vector3i size = end - begin;
            const auto count = static_cast<std::size_t>(size.prod());
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, static_cast<Eigen::Index>(count));
            std::vector<std::size_t> targets(count);
            std::size_t local = 0;
            for (int i = begin[0]; i < end[0]; ++i)
                for (int j = begin[1]; j < end[1]; ++j)
                    for (int k = begin[2]; k < end[2]; ++k, ++local) {
                        points.col(static_cast<Eigen::Index>(local)) =
                            grid.origin + grid.steps * Eigen::Vector3d(i, j, k);
                        targets[local] = (static_cast<std::size_t>(i) *
                                              static_cast<std::size_t>(grid.counts[1]) +
                                          static_cast<std::size_t>(j)) *
                                             static_cast<std::size_t>(grid.counts[2]) +
                                         static_cast<std::size_t>(k);
                    }
            const Eigen::Vector3d center =
                grid.origin +
                grid.steps * (0.5 * (begin + end - Eigen::Vector3i::Ones()).cast<double>());
            std::vector<double> values(count);
            potential.evaluate_cluster(points, center, plan, values);
            for (std::size_t p = 0; p < count; ++p)
                output[targets[p]] = values[p];
            report.advance(count);
        } catch (...) {
            worker_errors.capture();
        }
    }
    worker_errors.rethrow_if_any();
}
} // namespace

Evaluator::Evaluator(const Molecule& molecule, EvaluationOptions options)
    : molecule_(&molecule), options_(options) {
    if (!options.block_size || options.threads < 0)
        throw std::invalid_argument("invalid evaluation options");
    detail::validate_screening_threshold(options.screening_threshold);
    if (!std::isfinite(options.potential_tolerance) || options.potential_tolerance < 0.0)
        throw std::invalid_argument("potential tolerance must be finite and nonnegative");
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

    // Progress is reported over all passes on the grid: the point path (if
    // any), then one pass per field evaluated on the grid path below.
    const auto on_grid_path = [&](const Property& key) {
        return key.kind == PropertyKind::Orbital ||
               (key.kind == PropertyKind::Density && backend != Backend::AO_MO);
    };
    const bool point_pass = uses_point_path(plan, backend);
    const bool potential_pass = plan.needs_electrostatic_potential();
    std::size_t passes = (point_pass ? 1 : 0) + (potential_pass ? 1 : 0);
    for (const auto& entry : result.fields_)
        passes += on_grid_path(entry.first) ? 1 : 0;
    std::size_t completed = 0; // points of the passes already finished
    ProgressCallback pass_progress;
    if (options_.progress)
        pass_progress = [&](std::size_t done, std::size_t) {
            options_.progress(completed + done, passes * n);
        };

    if (point_pass) {
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
                            result.fields_, pass_progress);
        completed += n;
    }

    if (potential_pass) {
        // The density matrix is that of the FCHK with either backend; with
        // AO-MO it reproduces the orbital density (checked at load time).
        const detail::ElectrostaticPotential potential(molecule_->evaluator_.basis(),
                                                       molecule_->density_matrix(Spin::Total),
                                                       molecule_->atoms());
        evaluate_electrostatic_potential(
            grid, n, options_.potential_tolerance, potential,
            result.fields_.at(density::total.electrostatic_potential()).components[0],
            pass_progress);
        completed += n;
    }

    // Grid path: orbitals always, densities only with AO-P.
    for (auto& [key, output] : result.fields_) {
        if (!on_grid_path(key)) continue;
        const Field field =
            key.kind == PropertyKind::Density
                ? molecule_
                      ->densities_blocked(grid, {spin_of(key.component)}, options_.block_size,
                                          options_.screening_threshold, pass_progress)
                      .at(spin_of(key.component))
                : molecule_
                      ->orbitals(grid, {key.orbital_index}, spin_of(key.channel), pass_progress)
                      .front();
        std::copy(field.values().begin(), field.values().end(), output.components[0].begin());
        completed += n;
    }
    return result;
}
} // namespace chem::kumo
