#include "chem/kumo/evaluator.hpp"
#include "evaluation_plan.hpp"
#include "automatic_backend.hpp"
#include "electrostatic_potential.hpp"
#include "ao_screening.hpp"
#include "grid_geometry.hpp"
#include "orbital_energy_fields.hpp"
#include "parallel_blocks.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

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

// The points of one evaluation: a validated grid, in Field order, or the
// columns of a matrix.
struct PointSet {
    const Grid* grid = nullptr;
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>* columns = nullptr;
    std::size_t count = 0;

    [[nodiscard]] Eigen::Matrix<double, 3, Eigen::Dynamic> block(std::size_t start,
                                                                 std::size_t size) const {
        if (grid) return detail::grid_points(*grid, start, size);
        return columns->middleCols(static_cast<Eigen::Index>(start),
                                   static_cast<Eigen::Index>(size));
    }
};

// Evaluation paths:
//  - point path: shared AO derivative evaluation on explicit point blocks,
//    used for all density derivatives and, with AO-MO or at arbitrary
//    points, also for densities;
//  - grid path (grids only): the validated Molecule::density /
//    Molecule::orbital kernels, used for AO-P densities and for orbitals;
//  - at arbitrary points, orbitals come from their own pass;
//  - the electrostatic potential and the orbital-energy fields come from
//    their own passes.
bool uses_point_path(const detail::ComponentPlan& component, Backend backend, bool grid) {
    return component.derivatives() || ((backend == Backend::AO_MO || !grid) && component.density);
}

bool uses_point_path(const detail::EvaluationPlan& plan, Backend backend, bool grid) {
    return std::any_of(plan.components.begin(), plan.components.end(),
                       [&](const auto& c) { return uses_point_path(c, backend, grid); });
}

// Fields evaluated on the grid path: orbitals always, densities with AO-P
// (`backend` is the resolved backend, AO_P or AO_MO).
bool on_grid_path(const Property& key, Backend backend, bool grid) {
    return grid && (key.kind == PropertyKind::Orbital ||
                    (key.kind == PropertyKind::Density && backend == Backend::AO_P));
}

bool is_orbital_energy_field(PropertyKind kind) {
    return kind == PropertyKind::AverageLocalIonizationEnergy ||
           kind == PropertyKind::LocalElectronAffinity ||
           kind == PropertyKind::LocalElectronAttachmentEnergy ||
           kind == PropertyKind::LocalElectronegativity || kind == PropertyKind::LocalHardness;
}

// Fields that the point path stores.
bool from_point_path(const Property& key, Backend backend, bool grid) {
    return !on_grid_path(key, backend, grid) && key.kind != PropertyKind::Orbital &&
           key.kind != PropertyKind::ElectrostaticPotential && !is_orbital_energy_field(key.kind);
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
    constexpr double pi = std::numbers::pi;
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
    constexpr double pi = std::numbers::pi;
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

// Localized orbital locator (Schmider and Becke 2000), LOL = t / (1 + t)
// with t = D0 / tau: tau summed over the channels given, D0 the Thomas-Fermi
// reference of ELF summed over the channels with positive density. As in
// Multiwfn, 0 with no density; where tau = 0 and D0 > 0 (t infinite) it is
// 1, the limit (Multiwfn returns 0 there).
double lol(std::initializer_list<const DensityDerivatives*> channels) {
    constexpr double pi = std::numbers::pi;
    static const double fermi = 0.3 * std::pow(6.0 * pi * pi, 2.0 / 3.0);
    double tau = 0.0, uniform = 0.0;
    for (const DensityDerivatives* d : channels) {
        tau += d->tau;
        if (d->value > 0.0) uniform += fermi * std::pow(d->value, 5.0 / 3.0);
    }
    if (!(uniform > 0.0)) return 0.0;
    if (!(tau > 0.0)) return 1.0;
    return uniform / (tau + uniform);
}

// Interaction region indicator (Lu and Chen 2021), |grad rho| / rho^1.1,
// capped at interaction_region_indicator_cap (also where rho <= 0).
double interaction_region_indicator(const DensityDerivatives& d) {
    if (!(d.value > 0.0)) return interaction_region_indicator_cap;
    const double value = d.gradient.norm() / std::pow(d.value, 1.1);
    return value < interaction_region_indicator_cap ? value : interaction_region_indicator_cap;
}

// Density overlap regions indicator (de Silva and Corminboeuf 2014),
// DORI = theta / (1 + theta) with
//   theta = |grad (|grad rho|^2 / rho^2)|^2 / (|grad rho|^2 / rho^2)^3,
// computed with q = grad rho / rho, independent of the scale of rho:
//   theta = 4 |(H / rho) q - |q|^2 q|^2 / |q|^6.
// At a critical point (q = 0) theta is infinite and DORI is 1, its limit;
// where rho <= 0 it is 0, its limit far from the molecule.
double density_overlap_regions_indicator(const DensityDerivatives& d) {
    if (!(d.value > 0.0)) return 0.0;
    const Eigen::Vector3d q = d.gradient / d.value;
    const double q2 = q.squaredNorm();
    if (!(q2 > 0.0)) return 1.0;
    const Eigen::Vector3d w = d.hessian * q / d.value - q2 * q;
    const double theta = 4.0 * w.squaredNorm() / (q2 * q2 * q2);
    if (!std::isfinite(theta)) return 1.0;
    return theta / (1.0 + theta);
}

using BlockDerivatives = std::map<Spin, std::vector<DensityDerivatives>>;

// Copy (or derive) one property of block results, starting at grid point
// `start`. Evaluator Hessians use row-major upper-triangle order: xx, xy, xz,
// yy, yz, zz.
void store_block(PropertyField& output, const Property& key, const BlockDerivatives& block,
                 std::size_t count, std::size_t start) {
    if ((key.kind == PropertyKind::Elf || key.kind == PropertyKind::Lol) &&
        key.component == DensityComponent::Total) {
        const auto& alpha = block.at(Spin::Alpha);
        const auto& beta = block.at(Spin::Beta);
        for (std::size_t local = 0; local < count; ++local)
            output.components[0][start + local] = key.kind == PropertyKind::Elf
                                                      ? elf({&alpha[local], &beta[local]})
                                                      : lol({&alpha[local], &beta[local]});
        return;
    }
    if (key.kind == PropertyKind::SpinPolarization) {
        const auto& total = block.at(Spin::Total);
        const auto& spin = block.at(Spin::Difference);
        for (std::size_t local = 0; local < count; ++local)
            output.components[0][start + local] =
                total[local].value > 0.0 ? spin[local].value / total[local].value : 0.0;
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
        case PropertyKind::Lol: output.components[0][p] = lol({&d}); break;
        case PropertyKind::InteractionRegionIndicator:
            output.components[0][p] = interaction_region_indicator(d);
            break;
        case PropertyKind::DensityOverlapRegionsIndicator:
            output.components[0][p] = density_overlap_regions_indicator(d);
            break;
        case PropertyKind::HamiltonianKineticEnergy:
            output.components[0][p] = d.tau - 0.25 * d.laplacian;
            break;
        case PropertyKind::PotentialEnergyDensity:
            output.components[0][p] = 0.25 * d.laplacian - 2.0 * d.tau;
            break;
        case PropertyKind::EnergyDensity:
            output.components[0][p] = 0.25 * d.laplacian - d.tau;
            break;
        case PropertyKind::ShannonEntropyDensity:
            output.components[0][p] = d.value > 0.0 ? -d.value * std::log(d.value) : 0.0;
            break;
        case PropertyKind::FisherInformationDensity:
            output.components[0][p] = d.value > 0.0 ? d.gradient.squaredNorm() / d.value : 0.0;
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
void evaluate_point_path(const PointSet& set, std::size_t block_size,
                         const detail::EvaluationPlan& plan, Backend backend,
                         const BlockEvaluator& evaluate_block,
                         std::map<Property, PropertyField>& fields,
                         const ProgressCallback& progress) {
    const bool grid = set.grid != nullptr;
    const std::size_t n = set.count;
    std::vector<Spin> needed;
    for (auto component : all_components)
        if (uses_point_path(plan_for(plan, component), backend, grid))
            needed.push_back(spin_of(component));

    const std::size_t nblocks =
        detail::parallel_block_count(n, block_size, "too many point blocks for evaluation");
    detail::Progress report(progress, n);
    detail::for_each_block<detail::Schedule::Dynamic>(nblocks, [&](std::size_t index) {
        const std::size_t start = static_cast<std::size_t>(index) * block_size;
        const std::size_t count = std::min(block_size, n - start);
        const auto points = set.block(start, count);
        const auto block = evaluate_block(points, needed);

        for (auto& [key, output] : fields)
            if (from_point_path(key, backend, grid)) store_block(output, key, block, count, start);
        // Serialized and monotonic; a throwing callback is captured
        // like any other block error.
        report.advance(count);
    });
}

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

    const auto ntiles = static_cast<std::size_t>(tiles[0]) * static_cast<std::size_t>(tiles[1]) *
                        static_cast<std::size_t>(tiles[2]);
    detail::Progress report(progress, n);
    detail::for_each_block<detail::Schedule::Dynamic>(ntiles, [&](std::size_t index) {
        const auto per_plane = static_cast<std::size_t>(tiles[1] * tiles[2]);
        const Eigen::Vector3i tile(static_cast<int>(index / per_plane),
                                   static_cast<int>((index / static_cast<std::size_t>(tiles[2])) %
                                                    static_cast<std::size_t>(tiles[1])),
                                   static_cast<int>(index % static_cast<std::size_t>(tiles[2])));
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
                    targets[local] = detail::grid_index(grid, i, j, k);
                }
        const Eigen::Vector3d center =
            grid.origin +
            grid.steps * (0.5 * (begin + end - Eigen::Vector3i::Ones()).cast<double>());
        std::vector<double> values(count);
        potential.evaluate_cluster(points, center, plan, values);
        for (std::size_t p = 0; p < count; ++p)
            output[targets[p]] = values[p];
        report.advance(count);
    });
}

// Electrostatic potential at arbitrary points: points are grouped by the
// cube of edge potential_cluster_edge, aligned with the axes at multiples of
// the edge, that contains them, and each group is evaluated around its cube
// centre with one far-field plan covering a cube. Groups depend only on the
// points, so the result does not depend on the number of threads.
void evaluate_electrostatic_potential(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points, double tolerance,
    const detail::ElectrostaticPotential& potential, std::vector<double>& output,
    const ProgressCallback& progress) {
    const std::size_t n = static_cast<std::size_t>(points.cols());
    using Cell = std::array<std::int64_t, 3>;
    std::vector<Cell> cells(n);
    for (std::size_t p = 0; p < n; ++p)
        for (int a = 0; a < 3; ++a) {
            const double index =
                std::floor(points(a, static_cast<Eigen::Index>(p)) / potential_cluster_edge);
            if (!(std::abs(index) < 1e15))
                throw std::invalid_argument("point too far for the electrostatic potential");
            cells[p][static_cast<std::size_t>(a)] = static_cast<std::int64_t>(index);
        }
    std::vector<std::size_t> order(n);
    for (std::size_t p = 0; p < n; ++p)
        order[p] = p;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return cells[a] < cells[b]; });
    std::vector<std::size_t> starts;
    for (std::size_t p = 0; p < n; ++p)
        if (p == 0 || cells[order[p]] != cells[order[p - 1]]) starts.push_back(p);
    starts.push_back(n);
    // Half the cube diagonal, slightly enlarged to cover rounding.
    const double radius = 0.5 * std::sqrt(3.0) * potential_cluster_edge;
    const auto plan = potential.far_field_plan(radius * (1.0 + 1e-10) + 1e-12, tolerance);

    detail::Progress report(progress, n);
    detail::for_each_block<detail::Schedule::Dynamic>(starts.size() - 1, [&](std::size_t g) {
        const std::size_t begin = starts[g], end = starts[g + 1];
        const auto count = static_cast<Eigen::Index>(end - begin);
        Eigen::Matrix<double, 3, Eigen::Dynamic> cluster(3, count);
        for (Eigen::Index p = 0; p < count; ++p)
            cluster.col(p) =
                points.col(static_cast<Eigen::Index>(order[begin + static_cast<std::size_t>(p)]));
        const Cell& cell = cells[order[begin]];
        const Eigen::Vector3d center((static_cast<double>(cell[0]) + 0.5) * potential_cluster_edge,
                                     (static_cast<double>(cell[1]) + 0.5) * potential_cluster_edge,
                                     (static_cast<double>(cell[2]) + 0.5) * potential_cluster_edge);
        std::vector<double> values(static_cast<std::size_t>(count));
        potential.evaluate_cluster(cluster, center, plan, values);
        for (std::size_t p = 0; p < values.size(); ++p)
            output[order[begin + p]] = values[p];
        report.advance(values.size());
    });
}

// ALIE, LEA, LEAE, local electronegativity and hardness from the orbital
// sums (see detail::OrbitalEnergyFields), 0 where a denominator vanishes.
void store_orbital_energy_block(std::map<Property, PropertyField>& fields,
                                const detail::OrbitalEnergySums& sums, std::size_t count,
                                std::size_t start) {
    for (auto& [key, output] : fields) {
        if (!is_orbital_energy_field(key.kind)) continue;
        for (std::size_t t = 0; t < count; ++t) {
            const double rho = sums.occupied_density[t];
            const double alie = rho > 0.0 ? sums.occupied_energy[t] / rho : 0.0;
            double lea = 0.0, attachment = 0.0;
            if (!sums.virtual_density.empty()) {
                const double vir = sums.virtual_density[t];
                lea = vir > 0.0 ? -sums.virtual_energy[t] / vir : 0.0;
                attachment = rho > 0.0 ? sums.attachment_energy[t] / rho : 0.0;
            }
            double value = 0.0;
            switch (key.kind) {
            case PropertyKind::AverageLocalIonizationEnergy: value = alie; break;
            case PropertyKind::LocalElectronAffinity: value = lea; break;
            case PropertyKind::LocalElectronAttachmentEnergy: value = attachment; break;
            case PropertyKind::LocalElectronegativity: value = 0.5 * (alie + lea); break;
            default: value = 0.5 * (alie - lea); break;
            }
            output.components[0][start + t] = value;
        }
    }
}
} // namespace

Evaluator::Evaluator(const Molecule& molecule, EvaluationOptions options)
    : molecule_(&molecule), options_(options) {
    if (options.block_size == 0 || options.block_size > max_block_size)
        throw std::invalid_argument("evaluation block_size must be in [1, " +
                                    std::to_string(max_block_size) + "]");
    if (options.threads < 0) throw std::invalid_argument("invalid evaluation options");
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
    return evaluate_at(&grid, nullptr, n, request);
}

Results
Evaluator::evaluate(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                    const Request& request) const {
    if (!points.allFinite()) throw std::invalid_argument("evaluation points must be finite");
    return evaluate_at(nullptr, &points, static_cast<std::size_t>(points.cols()), request);
}

Results
Evaluator::evaluate_at(const Grid* grid,
                       const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>* points,
                       std::size_t n, const Request& request) const {
    const PointSet set{grid, points, n};
    const bool on_grid = grid != nullptr;
    // Plan the requested properties before resolving the electronic backend.
    const auto plan = detail::EvaluationPlan::from(request);
    const auto& info = molecule_->electronic_structure_info();
    const Backend backend = detail::choose_backend(options_.backend, info.ao_mo_eligible(), plan);
    Results result;
    result.backend_used_ = backend;
    for (auto key : request.properties()) {
        result.fields_.try_emplace(
            key, PropertyField{n, std::vector<std::vector<double>>(channels(key.kind),
                                                                   std::vector<double>(n))});
    }
    // Check the orbitals (and the orbital indices) before any pass.
    std::optional<detail::OrbitalEnergyFields> orbital_energies;
    if (plan.orbital_energies) {
        const char* what = plan.virtual_orbitals ? "the local electron affinity and related fields"
                                                 : "the average local ionization energy";
        orbital_energies.emplace(
            *molecule_, molecule_->evaluator_, molecule_->mo_coefficients(Spin::Alpha),
            molecule_->mo_coefficients(Spin::Beta), plan.virtual_orbitals, what);
    }
    std::vector<Property> orbitals; // at arbitrary points: their own pass
    for (const auto& entry : result.fields_)
        if (entry.first.kind == PropertyKind::Orbital) {
            const auto& coefficients = molecule_->mo_coefficients(spin_of(entry.first.channel));
            if (coefficients.cols() == 0) throw std::runtime_error("FCHK has no MO coefficients");
            if (entry.first.orbital_index > static_cast<std::size_t>(coefficients.cols()))
                throw std::out_of_range("MO index outside available coefficient range");
            if (!on_grid) orbitals.push_back(entry.first);
        }

    // Progress is reported over all passes: the point path (if any), the
    // potential, the orbital-energy fields and the orbitals at arbitrary
    // points, then one pass per field evaluated on the grid path.
    const bool point_pass = uses_point_path(plan, backend, on_grid);
    const bool potential_pass = plan.needs_electrostatic_potential();
    std::size_t passes = (point_pass ? 1 : 0) + (potential_pass ? 1 : 0) +
                         (orbital_energies ? 1 : 0) + (orbitals.empty() ? 0 : 1);
    for (const auto& entry : result.fields_)
        passes += on_grid_path(entry.first, backend, on_grid) ? 1 : 0;
    std::size_t completed = 0; // points of the passes already finished
    ProgressCallback pass_progress;
    if (options_.progress)
        pass_progress = [&](std::size_t done, std::size_t) {
            options_.progress(completed + done, passes * n);
        };

    if (point_pass) {
        // The AO-MO and AO-P block kernels are private to Molecule; Evaluator
        // is a friend, so the call is wrapped here and passed to the helper.
        // At arbitrary points densities come from this path, and the AO-P
        // kernel starts at first derivatives.
        const int ao_order = on_grid ? plan.maximum_ao_order : std::max(1, plan.maximum_ao_order);
        const auto evaluate_block = [&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& block,
                                        const std::vector<Spin>& needed) {
            return backend == Backend::AO_MO
                       ? molecule_->density_derivatives_mo(block, needed, options_.block_size,
                                                           std::max(0, ao_order),
                                                           options_.screening_threshold)
                       : molecule_->density_derivatives_shared(block, needed, options_.block_size,
                                                               ao_order,
                                                               options_.screening_threshold);
        };
        evaluate_point_path(set, options_.block_size, plan, backend, evaluate_block, result.fields_,
                            pass_progress);
        completed += n;
    }

    if (potential_pass) {
        // The density matrix is that of the FCHK with either backend; with
        // AO-MO it reproduces the orbital density (checked at load time).
        const detail::ElectrostaticPotential potential(molecule_->evaluator_.basis(),
                                                       molecule_->density_matrix(Spin::Total),
                                                       molecule_->atoms());
        auto& output = result.fields_.at(density::total.electrostatic_potential()).components[0];
        if (on_grid)
            evaluate_electrostatic_potential(*grid, n, options_.potential_tolerance, potential,
                                             output, pass_progress);
        else
            evaluate_electrostatic_potential(*points, options_.potential_tolerance, potential,
                                             output, pass_progress);
        completed += n;
    }

    // Blocks of points for the orbital passes: each block writes its own
    // range, so the result does not depend on the number of threads.
    const auto for_each_point_block = [&](const auto& body) {
        const std::size_t size = options_.block_size;
        const std::size_t nblocks =
            detail::parallel_block_count(n, size, "too many point blocks for evaluation");
        detail::Progress report(pass_progress, n);
        detail::for_each_block<detail::Schedule::Dynamic>(nblocks, [&](std::size_t index) {
            const std::size_t start = index * size;
            const std::size_t count = std::min(size, n - start);
            body(set.block(start, count), start, count);
            report.advance(count);
        });
        completed += n;
    };

    if (orbital_energies)
        for_each_point_block([&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& block,
                                 std::size_t start, std::size_t count) {
            store_orbital_energy_block(result.fields_, orbital_energies->evaluate(block), count,
                                       start);
        });

    if (!orbitals.empty())
        for_each_point_block([&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& block,
                                 std::size_t start, std::size_t count) {
            const auto& basis = molecule_->evaluator_;
            Eigen::MatrixXd phi(static_cast<Eigen::Index>(basis.basis_size()),
                                static_cast<Eigen::Index>(count));
            basis.evaluate_batch(block, phi);
            for (const Property& key : orbitals) {
                const auto& coefficients = molecule_->mo_coefficients(spin_of(key.channel));
                const Eigen::VectorXd psi =
                    phi.transpose() *
                    coefficients.col(static_cast<Eigen::Index>(key.orbital_index - 1));
                auto& output = result.fields_.at(key).components[0];
                for (std::size_t t = 0; t < count; ++t)
                    output[start + t] = psi[static_cast<Eigen::Index>(t)];
            }
        });

    // Grid path: orbitals and AO-P densities on grids.
    for (auto& [key, output] : result.fields_) {
        if (!on_grid_path(key, backend, on_grid)) continue;
        const Field field =
            key.kind == PropertyKind::Density
                ? molecule_
                      ->densities_blocked(*grid, {spin_of(key.component)}, options_.block_size,
                                          options_.screening_threshold, pass_progress)
                      .at(spin_of(key.component))
                : molecule_
                      ->orbitals(*grid, {key.orbital_index}, spin_of(key.channel), pass_progress)
                      .front();
        std::copy(field.values().begin(), field.values().end(), output.components[0].begin());
        completed += n;
    }
    return result;
}
} // namespace chem::kumo
