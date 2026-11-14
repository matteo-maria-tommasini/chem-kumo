#include "chem/kumo/evaluator.hpp"

#include "ao_screening.hpp"
#include "automatic_backend.hpp"
#include "electrostatic_potential.hpp"
#include "evaluation_plan.hpp"
#include "grid_geometry.hpp"
#include "orbital_energy_fields.hpp"
#include "parallel_blocks.hpp"
#include "potential_clusters.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
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
std::size_t channels(PropertyKind k) {
    switch (k) {
    case PropertyKind::Gradient:
    case PropertyKind::OrbitalGradient: return 3;
    case PropertyKind::Hessian:
    case PropertyKind::OrbitalHessian: return 6;
    default: return 1;
    }
}

// Properties of one molecular orbital: its amplitude and derivatives.
bool is_orbital_property(PropertyKind k) {
    return k == PropertyKind::Orbital || k == PropertyKind::OrbitalGradient ||
           k == PropertyKind::OrbitalHessian;
}

constexpr DensityComponent all_components[] = {DensityComponent::Total, DensityComponent::Alpha,
                                               DensityComponent::Beta,
                                               DensityComponent::Difference};

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

// Evaluation paths (AO = atomic orbital; the AO-P backend contracts AO values
// with the density matrix, AO-MO with the occupied molecular orbitals):
//  - point path: one shared evaluation of AO values and derivatives on
//    explicit blocks of points, used for all density derivatives and, with
//    AO-MO or at arbitrary points, also for densities;
//  - grid path (grids only): the validated Molecule::densities_blocked and
//    Molecule::orbitals kernels, used for AO-P densities and for orbitals;
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

// Fields evaluated on the grid path: orbital amplitudes (not their
// derivatives) always, densities with AO-P
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

// Whether this field is computed and stored by the point path.
bool from_point_path(const Property& key, Backend backend, bool grid) {
    return !on_grid_path(key, backend, grid) && !is_orbital_property(key.kind) &&
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
double elf(std::initializer_list<const detail::DensityDerivatives*> channels) {
    constexpr double pi = std::numbers::pi;
    static const double fermi = 0.3 * std::pow(6.0 * pi * pi, 2.0 / 3.0);
    double pauli = 0.0;
    double uniform = 0.0;
    for (const detail::DensityDerivatives* d : channels) {
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
double reduced_gradient(const detail::DensityDerivatives& d) {
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
double sign_lambda2_density(const detail::DensityDerivatives& d) {
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(d.hessian, Eigen::EigenvaluesOnly);
    return solver.eigenvalues()[1] < 0.0 ? -d.value : d.value;
}

// Localized orbital locator (Schmider and Becke 2000), LOL = t / (1 + t)
// with t = D0 / tau: tau summed over the channels given, D0 the Thomas-Fermi
// reference of ELF summed over the channels with positive density. As in
// Multiwfn, 0 with no density; where tau = 0 and D0 > 0 (t infinite) it is
// 1, the limit (Multiwfn returns 0 there).
double lol(std::initializer_list<const detail::DensityDerivatives*> channels) {
    constexpr double pi = std::numbers::pi;
    static const double fermi = 0.3 * std::pow(6.0 * pi * pi, 2.0 / 3.0);
    double tau = 0.0, uniform = 0.0;
    for (const detail::DensityDerivatives* d : channels) {
        tau += d->tau;
        if (d->value > 0.0) uniform += fermi * std::pow(d->value, 5.0 / 3.0);
    }
    if (!(uniform > 0.0)) return 0.0;
    if (!(tau > 0.0)) return 1.0;
    return uniform / (tau + uniform);
}

// Interaction region indicator (Lu and Chen 2021), |grad rho| / rho^1.1,
// capped at interaction_region_indicator_cap (also where rho <= 0).
double interaction_region_indicator(const detail::DensityDerivatives& d) {
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
double density_overlap_regions_indicator(const detail::DensityDerivatives& d) {
    if (!(d.value > 0.0)) return 0.0;
    const Eigen::Vector3d q = d.gradient / d.value;
    const double q2 = q.squaredNorm();
    if (!(q2 > 0.0)) return 1.0;
    const Eigen::Vector3d w = d.hessian * q / d.value - q2 * q;
    const double theta = 4.0 * w.squaredNorm() / (q2 * q2 * q2);
    if (!std::isfinite(theta)) return 1.0;
    return theta / (1.0 + theta);
}

using BlockDerivatives = std::map<DensityComponent, std::vector<detail::DensityDerivatives>>;

// Copy (or derive) one property from the density derivatives of a block of
// `count` points into `output`, starting at point index `start`. Evaluator
// Hessians use row-major upper-triangle order: xx, xy, xz, yy, yz, zz.
void store_block(PropertyField& output, const Property& key, const BlockDerivatives& block,
                 std::size_t count, std::size_t start) {
    if ((key.kind == PropertyKind::Elf || key.kind == PropertyKind::Lol) &&
        key.component == DensityComponent::Total) {
        const auto& alpha = block.at(DensityComponent::Alpha);
        const auto& beta = block.at(DensityComponent::Beta);
        for (std::size_t local = 0; local < count; ++local)
            output.components[0][start + local] = key.kind == PropertyKind::Elf
                                                      ? elf({&alpha[local], &beta[local]})
                                                      : lol({&alpha[local], &beta[local]});
        return;
    }
    if (key.kind == PropertyKind::SpinPolarization) {
        const auto& total = block.at(DensityComponent::Total);
        const auto& spin = block.at(DensityComponent::Difference);
        for (std::size_t local = 0; local < count; ++local)
            output.components[0][start + local] =
                total[local].value > 0.0 ? spin[local].value / total[local].value : 0.0;
        return;
    }
    const auto& derivatives = block.at(key.component);
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

// Point path. Points are processed in blocks of at most block_size, so
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
    std::vector<DensityComponent> needed;
    for (auto component : all_components)
        if (uses_point_path(plan_for(plan, component), backend, grid)) needed.push_back(component);

    detail::Progress report(progress, n);
    detail::for_each_point_block(
        n, block_size, "too many point blocks for evaluation", [&](detail::PointBlock range) {
            const auto block = evaluate_block(set.block(range.start, range.size), needed);
            for (auto& [key, output] : fields)
                if (from_point_path(key, backend, grid))
                    store_block(output, key, block, range.size, range.start);
            // Progress reports are serialized and never decrease; an
            // exception thrown by the callback is captured like any other
            // block error.
            report.advance(range.size);
        });
}

// Electrostatic potential pass on a grid. The grid is split into tiles of up
// to cluster_side^3 points (boxes of grid indices), each evaluated as one
// target cluster around the centre of its box, with one far-field plan whose
// radius covers the largest tile (see 'Target clusters' in
// docs/mathematics/chem-kumo-mathematics.tex). Cluster
// boundaries depend only on the grid, and each value only on its point, its
// cluster and the plan, so the result does not depend on the number of
// threads or on the block size.
constexpr int cluster_side = 8;

void evaluate_electrostatic_potential(const Grid& grid, std::size_t n, double tolerance,
                                      const detail::ElectrostaticPotential& potential,
                                      std::vector<double>& output,
                                      const ProgressCallback& progress) {
    Eigen::Vector3i tiles;
    Eigen::Vector3d edges; // of the largest tile, in grid steps
    for (int a = 0; a < 3; ++a) {
        edges[a] = std::min(cluster_side, grid.counts[a]) - 1;
        tiles[a] = (grid.counts[a] + cluster_side - 1) / cluster_side;
    }
    const auto plan =
        detail::cluster_plan(potential, detail::grid_box_radius(grid.steps, edges), tolerance);
    const auto ntiles = static_cast<std::size_t>(tiles[0]) * static_cast<std::size_t>(tiles[1]) *
                        static_cast<std::size_t>(tiles[2]);
    const auto gather = [&](std::size_t index, Eigen::Matrix<double, 3, Eigen::Dynamic>& points,
                            std::vector<std::size_t>& targets) {
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
        const auto count = static_cast<std::size_t>((end - begin).prod());
        points.resize(3, static_cast<Eigen::Index>(count));
        targets.resize(count);
        std::size_t local = 0;
        for (int i = begin[0]; i < end[0]; ++i)
            for (int j = begin[1]; j < end[1]; ++j)
                for (int k = begin[2]; k < end[2]; ++k, ++local) {
                    points.col(static_cast<Eigen::Index>(local)) =
                        grid.origin + grid.steps * Eigen::Vector3d(i, j, k);
                    targets[local] = detail::grid_index(grid, i, j, k);
                }
        return Eigen::Vector3d(grid.origin +
                               grid.steps *
                                   (0.5 * (begin + end - Eigen::Vector3i::Ones()).cast<double>()));
    };
    detail::evaluate_in_clusters(potential, plan, ntiles, gather, output, n, progress);
}

// Electrostatic potential at arbitrary points: points are grouped by the
// cube of edge potential_cluster_edge, aligned with the axes at multiples of
// the edge, that contains them, and each group is evaluated around its cube
// centre with one far-field plan covering a cube.
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
    const auto clusters = detail::cluster_by(n, [&](std::size_t p) { return cells[p]; });
    // Half the cube diagonal.
    const auto plan =
        detail::cluster_plan(potential, 0.5 * std::sqrt(3.0) * potential_cluster_edge, tolerance);
    const auto gather = [&](std::size_t c, Eigen::Matrix<double, 3, Eigen::Dynamic>& cluster,
                            std::vector<std::size_t>& targets) {
        detail::gather_cluster(
            clusters, c, [&](std::size_t p) { return points.col(static_cast<Eigen::Index>(p)); },
            cluster, targets);
        const Cell& cell = cells[targets.front()];
        return Eigen::Vector3d((static_cast<double>(cell[0]) + 0.5) * potential_cluster_edge,
                               (static_cast<double>(cell[1]) + 0.5) * potential_cluster_edge,
                               (static_cast<double>(cell[2]) + 0.5) * potential_cluster_edge);
    };
    detail::evaluate_in_clusters(potential, plan, clusters.size(), gather, output, n, progress);
}

// Average local ionization energy (ALIE), local electron affinity (LEA),
// local electron attachment energy (LEAE), local electronegativity and
// hardness from the orbital sums (see detail::OrbitalEnergySums), 0 where a
// denominator vanishes.
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

// One evaluation: the points, the plan of the request, the backend, the
// output fields and the progress over all passes. The passes are
//  - the point path (if any): density derivatives and, with AO-MO or at
//    arbitrary points, densities;
//  - the electrostatic potential;
//  - the orbital-energy fields;
//  - the orbitals at arbitrary points;
//  - the grid path: one pass per AO-P density or orbital on a grid.
// Progress is reported over all of them, in this order.
class Evaluator::Run {
public:
    Run(const Evaluator& evaluator, const Grid* grid,
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>* points, std::size_t n,
        const Request& request)
        : molecule_(*evaluator.molecule_), options_(evaluator.options_), set_{grid, points, n},
          plan_(detail::EvaluationPlan::from(request)) {
        // The backend choice depends on the plan (built above): AUTO picks
        // AO_MO only if some property needs a density or its derivatives.
        backend_ = detail::choose_backend(
            options_.backend, molecule_.electronic_structure_info().ao_mo_eligible(), plan_);
        result_.backend_used_ = backend_;
        for (auto key : request.properties())
            result_.fields_.try_emplace(
                key, PropertyField{n, std::vector<std::vector<double>>(channels(key.kind),
                                                                       std::vector<double>(n))});
    }

    Run(const Run&) = delete;
    Run& operator=(const Run&) = delete;

    Results execute() {
        // Check the orbitals (and the orbital indices) before any pass.
        std::optional<detail::OrbitalEnergyFields> orbital_energies;
        if (plan_.orbital_energies) {
            const char* what = plan_.virtual_orbitals
                                   ? "the local electron affinity and related fields"
                                   : "the average local ionization energy";
            orbital_energies.emplace(
                molecule_, molecule_.evaluator_, molecule_.mo_coefficients(SpinChannel::Alpha),
                molecule_.mo_coefficients(SpinChannel::Beta), plan_.virtual_orbitals, what);
        }
        const auto orbitals = checked_orbitals();

        const bool point_path = uses_point_path(plan_, backend_, on_grid());
        const bool potential = plan_.needs_electrostatic_potential();
        passes_ = (point_path ? 1 : 0) + (potential ? 1 : 0) + (orbital_energies ? 1 : 0) +
                  (orbitals.empty() ? 0 : 1);
        passes_ += grid_path_passes();
        if (options_.progress)
            pass_progress_ = [this](std::size_t done, std::size_t) {
                options_.progress(completed_ + done, passes_ * set_.count);
            };

        if (point_path) point_pass();
        if (potential) potential_pass();
        if (orbital_energies) orbital_energy_pass(*orbital_energies);
        if (!orbitals.empty()) orbital_pass(orbitals);
        grid_path();
        return std::move(result_);
    }

private:
    [[nodiscard]] bool on_grid() const { return set_.grid != nullptr; }

    // Check every requested orbital property; return those evaluated by the
    // orbital pass: all of them at arbitrary points, the derivatives on grids
    // (orbital amplitudes on grids take the grid path).
    [[nodiscard]] std::vector<Property> checked_orbitals() const {
        std::vector<Property> orbitals;
        for (const auto& entry : result_.fields_)
            if (is_orbital_property(entry.first.kind)) {
                const auto& coefficients = molecule_.mo_coefficients(entry.first.channel);
                if (coefficients.cols() == 0)
                    throw std::runtime_error("FCHK has no MO coefficients");
                if (entry.first.orbital_index > static_cast<std::size_t>(coefficients.cols()))
                    throw std::out_of_range("MO index outside available coefficient range");
                if (!on_grid_path(entry.first, backend_, on_grid()))
                    orbitals.push_back(entry.first);
            }
        return orbitals;
    }

    void point_pass() {
        // At arbitrary points densities come from this path too, and the AO-P
        // kernel always computes at least first derivatives, so ask for
        // order 1 or more.
        const int ao_order =
            on_grid() ? plan_.maximum_ao_order : std::max(1, plan_.maximum_ao_order);
        const auto evaluate_block = [&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& block,
                                        const std::vector<DensityComponent>& needed) {
            return backend_ == Backend::AO_MO
                       ? molecule_.density_derivatives_mo(block, needed, options_.block_size,
                                                          std::max(0, ao_order),
                                                          options_.screening_threshold)
                       : molecule_.density_derivatives_shared(block, needed, options_.block_size,
                                                              ao_order,
                                                              options_.screening_threshold);
        };
        evaluate_point_path(set_, options_.block_size, plan_, backend_, evaluate_block,
                            result_.fields_, pass_progress_);
        completed_ += set_.count;
    }

    void potential_pass() {
        // The potential always uses the FCHK total density matrix, whatever
        // the backend; when AO-MO is in use, that matrix was checked at load
        // time to agree with the occupied orbitals.
        const detail::ElectrostaticPotential potential(
            molecule_.evaluator_.basis(), molecule_.density_matrix(DensityComponent::Total),
            molecule_.atoms());
        auto& output = result_.fields_.at(density::total.electrostatic_potential()).components[0];
        if (on_grid())
            evaluate_electrostatic_potential(*set_.grid, set_.count, options_.potential_tolerance,
                                             potential, output, pass_progress_);
        else
            evaluate_electrostatic_potential(*set_.columns, options_.potential_tolerance, potential,
                                             output, pass_progress_);
        completed_ += set_.count;
    }

    void orbital_energy_pass(const detail::OrbitalEnergyFields& fields) {
        for_each_point_block(
            [&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& block, detail::PointBlock range) {
                store_orbital_energy_block(result_.fields_, fields.evaluate(block), range.size,
                                           range.start);
            });
    }

    // Orbital amplitudes and their derivatives on blocks of points. The
    // amplitudes come from the AO values; the derivatives from the AO
    // derivative panels (value, x, y, z, then the Hessian xx, xy, xz, yy, yz,
    // zz), computed only if a derivative is requested.
    void orbital_pass(const std::vector<Property>& orbitals) {
        const auto& basis = molecule_.evaluator_;
        const auto n = static_cast<Eigen::Index>(basis.basis_size());
        int derivative_order = 0;
        for (const Property& key : orbitals)
            if (key.kind == PropertyKind::OrbitalGradient)
                derivative_order = std::max(derivative_order, 1);
            else if (key.kind == PropertyKind::OrbitalHessian)
                derivative_order = 2;
        for_each_point_block(
            [&](const Eigen::Matrix<double, 3, Eigen::Dynamic>& block, detail::PointBlock range) {
                const auto b = static_cast<Eigen::Index>(range.size);
                Eigen::MatrixXd phi(n, b);
                basis.evaluate_batch(block, phi);
                Eigen::MatrixXd panels;
                if (derivative_order > 0) {
                    panels.resize(n, BasisEvaluator::derivative_panel_count(derivative_order) * b);
                    basis.evaluate_derivatives_batch(block, derivative_order, panels);
                }
                for (const Property& key : orbitals) {
                    const auto c = molecule_.mo_coefficients(key.channel)
                                       .col(static_cast<Eigen::Index>(key.orbital_index - 1));
                    auto& output = result_.fields_.at(key);
                    // Component `component` of the result from AO panel `panel`.
                    const auto store = [&](const Eigen::Ref<const Eigen::MatrixXd>& ao,
                                           std::size_t component) {
                        const Eigen::VectorXd psi = ao.transpose() * c;
                        for (std::size_t t = 0; t < range.size; ++t)
                            output.components[component][range.start + t] =
                                psi[static_cast<Eigen::Index>(t)];
                    };
                    if (key.kind == PropertyKind::Orbital)
                        store(phi, 0);
                    else if (key.kind == PropertyKind::OrbitalGradient)
                        for (std::size_t a = 0; a < 3; ++a)
                            store(panels.middleCols((1 + static_cast<Eigen::Index>(a)) * b, b), a);
                    else
                        for (std::size_t h = 0; h < 6; ++h)
                            store(panels.middleCols((4 + static_cast<Eigen::Index>(h)) * b, b), h);
                }
            });
    }

    // The grid-path fields of one kind: AO-P densities (one pass each) or
    // the orbital amplitudes of one spin channel (one pass for all of them,
    // sharing the AO values).
    [[nodiscard]] std::vector<Property> grid_path_fields(PropertyKind kind,
                                                         SpinChannel channel) const {
        std::vector<Property> keys;
        for (const auto& entry : result_.fields_)
            if (on_grid_path(entry.first, backend_, on_grid()) && entry.first.kind == kind &&
                (kind != PropertyKind::Orbital || entry.first.channel == channel))
                keys.push_back(entry.first);
        return keys;
    }

    [[nodiscard]] std::size_t grid_path_passes() const {
        std::size_t passes = grid_path_fields(PropertyKind::Density, SpinChannel::Alpha).size();
        for (const SpinChannel channel : {SpinChannel::Alpha, SpinChannel::Beta})
            passes += grid_path_fields(PropertyKind::Orbital, channel).empty() ? 0 : 1;
        return passes;
    }

    // AO-P densities and orbital amplitudes on grids, with the validated
    // Molecule kernels.
    void grid_path() {
        const auto copy = [&](const Property& key, const Field& field) {
            std::copy(field.values().begin(), field.values().end(),
                      result_.fields_.at(key).components[0].begin());
        };
        for (const Property& key : grid_path_fields(PropertyKind::Density, SpinChannel::Alpha)) {
            copy(key, molecule_
                          .densities_blocked(*set_.grid, {key.component}, options_.block_size,
                                             options_.screening_threshold, pass_progress_)
                          .at(key.component));
            completed_ += set_.count;
        }
        for (const SpinChannel channel : {SpinChannel::Alpha, SpinChannel::Beta}) {
            const auto keys = grid_path_fields(PropertyKind::Orbital, channel);
            if (keys.empty()) continue;
            std::vector<std::size_t> indices;
            for (const Property& key : keys)
                indices.push_back(key.orbital_index);
            const auto fields = molecule_.orbitals(*set_.grid, indices, channel, pass_progress_);
            for (std::size_t i = 0; i < keys.size(); ++i)
                copy(keys[i], fields[i]);
            completed_ += set_.count;
        }
    }

    // One pass over blocks of points, calling body(points, range) for each
    // block. Each block writes only its own range of points, so the result
    // does not depend on the number of threads.
    template <typename Body> void for_each_point_block(const Body& body) {
        detail::Progress report(pass_progress_, set_.count);
        detail::for_each_point_block(set_.count, options_.block_size,
                                     "too many point blocks for evaluation",
                                     [&](detail::PointBlock range) {
                                         body(set_.block(range.start, range.size), range);
                                         report.advance(range.size);
                                     });
        completed_ += set_.count;
    }

    const Molecule& molecule_;
    const EvaluationOptions& options_;
    const PointSet set_;
    const detail::EvaluationPlan plan_;
    Backend backend_ = Backend::AO_P;
    Results result_;
    std::size_t passes_ = 0;
    std::size_t completed_ = 0; // points of the passes already finished
    ProgressCallback pass_progress_;
};

Results Evaluator::evaluate(const Grid& grid, const Request& request) const {
    // Validate the grid up front: the derivative paths evaluate on explicit
    // points and would otherwise accept a degenerate geometry.
    const std::size_t n = detail::validated_point_count(grid);
    return Run(*this, &grid, nullptr, n, request).execute();
}

Results
Evaluator::evaluate(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                    const Request& request) const {
    if (!points.allFinite()) throw std::invalid_argument("evaluation points must be finite");
    return Run(*this, nullptr, &points, static_cast<std::size_t>(points.cols()), request).execute();
}
} // namespace chem::kumo
