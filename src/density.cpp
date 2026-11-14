// Electron densities and their derivatives, from atomic-orbital (AO) values:
// contractions with the FCHK AO density matrices P (AO-P), and the validated
// contraction with the occupied molecular orbitals (AO-MO).

#include "chem/kumo/density.hpp"

#include "ao_screening.hpp"
#include "density_matrix.hpp"
#include "derivative_layout.hpp"
#include "grid_geometry.hpp"
#include "parallel_blocks.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <utility>

namespace chem::kumo {
namespace {
// Density derivatives at point j of a block of b points. `ao` holds AO panels
// of b columns each: values, then the x, y, z first derivatives, then, for
// ao_order 2, the six second-derivative panels in ao_hessian_column order.
// `product` holds P times the first four panels (values and gradients).
// `out` must start zeroed; Hessians are left zero for ao_order 1.
void contract_derivatives(const Eigen::Ref<const Eigen::MatrixXd>& ao,
                          const Eigen::Ref<const Eigen::MatrixXd>& product, Eigen::Index b,
                          Eigen::Index j, int ao_order, DensityDerivatives& out) {
    const auto y = product.col(j);
    out.value = ao.col(j).dot(y);
    for (int a = 0; a < 3; ++a) {
        const auto ga = ao.col((a + 1) * b + j);
        out.gradient[a] = 2.0 * ga.dot(y);
        out.tau += 0.5 * ga.dot(product.col((a + 1) * b + j));
        if (ao_order == 2)
            for (int c = a; c < 3; ++c) {
                const auto hac = ao.col((4 + detail::ao_hessian_column[a][c]) * b + j);
                const double v = 2.0 * (hac.dot(y) + ga.dot(product.col((c + 1) * b + j)));
                out.hessian(a, c) = out.hessian(c, a) = v;
            }
    }
    out.laplacian = out.hessian.trace();
}

// scale * (a + sign * b), quantity by quantity (all are linear in P). With
// scale 1/2 this gives alpha = (total + spin) / 2 and beta = (total - spin) / 2
// from the total and spin contractions; with scale 1, total = alpha + beta and
// spin = alpha - beta. The Laplacian is the trace of the combined Hessian.
DensityDerivatives combine(const DensityDerivatives& a, const DensityDerivatives& b, double sign,
                           double scale) {
    DensityDerivatives d;
    d.value = scale * (a.value + sign * b.value);
    d.gradient = scale * (a.gradient + sign * b.gradient);
    d.hessian = scale * (a.hessian + sign * b.hessian);
    d.laplacian = d.hessian.trace();
    d.tau = scale * (a.tau + sign * b.tau);
    return d;
}

// Store `value` at `slot` of the component `spin`, if it was requested.
void store(std::map<DensityComponent, std::vector<DensityDerivatives>>& result,
           DensityComponent spin, std::size_t slot, const DensityDerivatives& value) {
    if (auto it = result.find(spin); it != result.end()) it->second[slot] = value;
}
} // namespace

DensityDerivatives Molecule::density_derivatives(const Eigen::Vector3d& point,
                                                 DensityComponent component) const {
    if (!point.allFinite()) throw std::invalid_argument("density evaluation point must be finite");
    // P is the (symmetric) AO density matrix of the requested component; for
    // alpha and beta it is built from the FCHK total and spin matrices as
    // (P_total +/- P_spin) / 2, so the same contraction serves all four.
    Eigen::MatrixXd storage;
    const Eigen::MatrixXd& p = detail::select_density_matrix(total_, spin_, component, storage);
    const auto ao = evaluator_.evaluate_derivatives(point);
    const Eigen::VectorXd y = p * ao.values;
    const Eigen::Matrix<double, Eigen::Dynamic, 3> pg = p * ao.gradient;
    DensityDerivatives result;
    result.value = ao.values.dot(y);
    result.gradient = 2.0 * ao.gradient.transpose() * y;
    for (int a = 0; a < 3; ++a) {
        for (int b = a; b < 3; ++b) {
            const double value = 2.0 * (ao.hessian.col(detail::ao_hessian_column[a][b]).dot(y) +
                                        ao.gradient.col(a).dot(pg.col(b)));
            result.hessian(a, b) = value;
            result.hessian(b, a) = value;
        }
        result.tau += 0.5 * ao.gradient.col(a).dot(pg.col(a));
    }
    result.laplacian = result.hessian.trace();
    return result;
}

std::vector<DensityDerivatives> Molecule::density_derivatives_blocked(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    DensityComponent component, std::size_t block_size, double screening_threshold) const {
    if (block_size == 0) throw std::invalid_argument("block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("density points must be finite");
    detail::validate_screening_threshold(screening_threshold);
    Eigen::MatrixXd storage;
    const Eigen::MatrixXd& p = detail::select_density_matrix(total_, spin_, component, storage);
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto count = static_cast<std::size_t>(points.cols());
    const auto capacity = static_cast<Eigen::Index>(std::min(block_size, count));
    std::vector<DensityDerivatives> result(count);

    struct Workspace {
        Eigen::MatrixXd ao; // value, x, y, z, then the AO Hessian (xx, yy, zz, xy, xz, yz)
        Eigen::MatrixXd product;
        Eigen::MatrixXd ao_kept, p_kept;
    };

    // Blocks write to disjoint output slots; each worker reuses its panels.
    detail::for_each_point_block_with_workspace(
        count, block_size, "too many derivative blocks for parallel evaluation",
        [&] {
            return Workspace{
                Eigen::MatrixXd(n, 10 * capacity), Eigen::MatrixXd(n, 4 * capacity), {}, {}};
        },
        [&](detail::PointBlock range, Workspace& w) {
            const auto b = static_cast<Eigen::Index>(range.size);
            const auto all = w.ao.leftCols(10 * b);
            evaluator_.evaluate_derivatives_batch(
                points.middleCols(static_cast<Eigen::Index>(range.start), b), 2, all);
            // Contract only the AOs that are significant somewhere in the block.
            const auto kept = detail::significant_aos(all, screening_threshold);
            const auto ao = detail::kept_rows(all, kept, w.ao_kept);
            const Eigen::MatrixXd& pk = detail::kept_submatrix(p, kept, w.p_kept);
            const auto m = static_cast<Eigen::Index>(kept.size());
            // The value and gradient panels are contiguous, so a single
            // matrix-matrix product applies P to all four, reading the
            // (potentially large) P only once.
            auto product = w.product.topLeftCorner(m, 4 * b);
            product.noalias() = pk * ao.leftCols(4 * b);
            for (Eigen::Index j = 0; j < b; ++j)
                contract_derivatives(ao, product, b, j, 2,
                                     result[range.start + static_cast<std::size_t>(j)]);
        });
    return result;
}

// One AO evaluation shared by all requested SCF spin components.
// Total and spin are contracted with the FCHK matrices; alpha and beta are
// linear combinations of these contractions, never reconstructed from MO
// occupations.
std::map<DensityComponent, std::vector<DensityDerivatives>> Molecule::density_derivatives_shared(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    const std::vector<DensityComponent>& components, std::size_t block_size, int ao_order,
    double screening_threshold) const {
    if (ao_order != 1 && ao_order != 2)
        throw std::invalid_argument("shared AO derivative order must be 1 or 2");
    if (block_size == 0) throw std::invalid_argument("block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("density points must be finite");
    detail::validate_screening_threshold(screening_threshold);
    const std::size_t count = static_cast<std::size_t>(points.cols());
    std::map<DensityComponent, std::vector<DensityDerivatives>> result;
    bool total_needed = false, spin_needed = false;
    for (DensityComponent component : components) {
        switch (component) {
        case DensityComponent::Total: total_needed = true; break;
        case DensityComponent::Difference: spin_needed = true; break;
        case DensityComponent::Alpha:
        case DensityComponent::Beta: total_needed = spin_needed = true; break;
        default: throw std::invalid_argument("invalid density spin component");
        }
        result.try_emplace(component, count);
    }
    if (result.empty() || count == 0) return result;
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto capacity = static_cast<Eigen::Index>(std::min(block_size, count));
    // AO panels: value, x, y, z and, for order 2, the AO Hessian.
    const int panels = BasisEvaluator::derivative_panel_count(ao_order);

    struct Workspace {
        Eigen::MatrixXd ao, total_product, spin_product;
        Eigen::MatrixXd ao_kept, total_kept, spin_kept;
    };

    const auto make_workspace = [&] {
        Workspace w;
        w.ao.resize(n, panels * capacity);
        if (total_needed) w.total_product.resize(n, 4 * capacity);
        if (spin_needed) w.spin_product.resize(n, 4 * capacity);
        return w;
    };
    detail::for_each_point_block_with_workspace(
        count, block_size, "too many derivative blocks for parallel evaluation", make_workspace,
        [&](detail::PointBlock range, Workspace& w) {
            const auto b = static_cast<Eigen::Index>(range.size);
            const auto all = w.ao.leftCols(panels * b);
            evaluator_.evaluate_derivatives_batch(
                points.middleCols(static_cast<Eigen::Index>(range.start), b), ao_order, all);
            // Contract only the AOs that are significant somewhere in the block.
            const auto kept = detail::significant_aos(all, screening_threshold);
            const auto ao = detail::kept_rows(all, kept, w.ao_kept);
            const auto m = static_cast<Eigen::Index>(kept.size());
            // P times the value and gradient panels (empty if not needed).
            const auto product_panel = [&](Eigen::MatrixXd& storage, bool needed) {
                return storage.topLeftCorner(needed ? m : 0, needed ? 4 * b : 0);
            };
            auto product_total = product_panel(w.total_product, total_needed);
            auto product_spin = product_panel(w.spin_product, spin_needed);
            if (total_needed)
                product_total.noalias() =
                    detail::kept_submatrix(total_, kept, w.total_kept) * ao.leftCols(4 * b);
            if (spin_needed)
                product_spin.noalias() =
                    detail::kept_submatrix(spin_, kept, w.spin_kept) * ao.leftCols(4 * b);
            // Contract the total and spin matrices once per point, then derive
            // alpha and beta from them as needed.
            for (Eigen::Index j = 0; j < b; ++j) {
                DensityDerivatives t, s;
                if (total_needed) contract_derivatives(ao, product_total, b, j, ao_order, t);
                if (spin_needed) contract_derivatives(ao, product_spin, b, j, ao_order, s);
                const auto slot = range.start + static_cast<std::size_t>(j);
                store(result, DensityComponent::Total, slot, t);
                store(result, DensityComponent::Difference, slot, s);
                if (result.contains(DensityComponent::Alpha))
                    store(result, DensityComponent::Alpha, slot, combine(t, s, 1.0, 0.5));
                if (result.contains(DensityComponent::Beta))
                    store(result, DensityComponent::Beta, slot, combine(t, s, -1.0, 0.5));
            }
        });
    return result;
}

// Validated occupied-MO contraction for electron-density derivatives.
// Evaluator selects this backend when the orbital representation is eligible;
// AO-P remains available as the general density-matrix reference.
std::map<DensityComponent, std::vector<DensityDerivatives>> Molecule::density_derivatives_mo(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    const std::vector<DensityComponent>& components, std::size_t block_size, int ao_order,
    double screening_threshold) const {
    if (!electronic_info_.ao_mo_eligible())
        throw std::invalid_argument("AO-MO backend is not eligible for this molecule");
    if (ao_order < 0 || ao_order > 2)
        throw std::invalid_argument("AO-MO derivative order must be 0, 1 or 2");
    if (block_size == 0 || block_size > max_block_size)
        throw std::invalid_argument("invalid AO-MO block size");
    if (!points.allFinite()) throw std::invalid_argument("nonfinite AO-MO points");
    detail::validate_screening_threshold(screening_threshold);
    const std::size_t count = static_cast<std::size_t>(points.cols());
    std::map<DensityComponent, std::vector<DensityDerivatives>> result;
    for (DensityComponent c : components) {
        switch (c) {
        case DensityComponent::Total:
        case DensityComponent::Alpha:
        case DensityComponent::Beta:
        case DensityComponent::Difference: break;
        default: throw std::invalid_argument("invalid AO-MO spin component");
        }
        result.try_emplace(c, count);
    }
    if (result.empty() || count == 0) return result;
    const bool need_a = result.contains(DensityComponent::Alpha) ||
                        result.contains(DensityComponent::Total) ||
                        result.contains(DensityComponent::Difference);
    const bool need_b = result.contains(DensityComponent::Beta) ||
                        result.contains(DensityComponent::Total) ||
                        result.contains(DensityComponent::Difference);
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto capacity = static_cast<Eigen::Index>(std::min(count, block_size));
    const auto& cb = unrestricted_ ? beta_mos_ : alpha_mos_;
    const Eigen::MatrixXd ca_occ = alpha_mos_.leftCols(nalpha_);
    const Eigen::MatrixXd cb_occ = cb.leftCols(nbeta_);
    // Closed shell: the beta occupied orbitals are the alpha ones, so the
    // beta projection and contraction would repeat the alpha ones exactly.
    const bool beta_is_alpha = need_a && need_b && !unrestricted_ && nalpha_ == nbeta_;
    // AO panels: value, then x, y, z and the AO Hessian (xx, yy, zz, xy, xz,
    // yz) as far as ao_order requires.
    const int panels = ao_order == 0 ? 1 : BasisEvaluator::derivative_panel_count(ao_order);
    using Projection = std::array<Eigen::MatrixXd, 10>;
    // Each spin channel is projected once per block onto its occupied MOs:
    // projected[a] = c^T * (AO panel a), with c the occupied MO coefficients,
    // both restricted to the AOs kept by screening.
    const auto project = [&](const Eigen::Ref<const Eigen::MatrixXd>& c,
                             const Eigen::Ref<const Eigen::MatrixXd>& panels_kept, Eigen::Index b,
                             Projection& projected) {
        if (c.cols() == 0) return;
        for (int a = 0; a < panels; ++a)
            projected[a].noalias() = c.transpose() * panels_kept.middleCols(a * b, b);
    };
    // Density derivatives of one spin channel at point j of the block.
    const auto contract = [&](const Projection& v, Eigen::Index j) {
        DensityDerivatives d;
        if (v[0].size() == 0) return d;
        d.value = v[0].col(j).squaredNorm();
        if (ao_order == 0) return d;
        for (int a = 0; a < 3; ++a) {
            d.gradient[a] = 2.0 * v[0].col(j).dot(v[a + 1].col(j));
            d.tau += 0.5 * v[a + 1].col(j).squaredNorm();
            if (ao_order == 2)
                for (int c = a; c < 3; ++c) {
                    const double value =
                        2.0 * (v[0].col(j).dot(v[4 + detail::ao_hessian_column[a][c]].col(j)) +
                               v[a + 1].col(j).dot(v[c + 1].col(j)));
                    d.hessian(a, c) = d.hessian(c, a) = value;
                }
        }
        d.laplacian = d.hessian.trace();
        return d;
    };

    struct Workspace {
        Eigen::MatrixXd ao, ao_kept, ca_kept, cb_kept;
    };

    detail::for_each_point_block_with_workspace(
        count, block_size, "too many AO-MO blocks",
        [&] { return Workspace{Eigen::MatrixXd(n, panels * capacity), {}, {}, {}}; },
        [&](detail::PointBlock range, Workspace& w) {
            const auto b = static_cast<Eigen::Index>(range.size);
            const auto block_points = points.middleCols(static_cast<Eigen::Index>(range.start), b);
            const auto all = w.ao.leftCols(panels * b);
            if (ao_order == 0)
                evaluator_.evaluate_batch(block_points, all);
            else
                evaluator_.evaluate_derivatives_batch(block_points, ao_order, all);
            // Screening: AOs negligible (value and derivatives) everywhere in
            // the block are dropped from the AO panels and the MO coefficients.
            const auto kept = detail::significant_aos(all, screening_threshold);
            const auto panels_kept = detail::kept_rows(all, kept, w.ao_kept);
            Projection va, vb;
            if (need_a) project(detail::kept_rows(ca_occ, kept, w.ca_kept), panels_kept, b, va);
            if (need_b && !beta_is_alpha)
                project(detail::kept_rows(cb_occ, kept, w.cb_kept), panels_kept, b, vb);
            for (Eigen::Index j = 0; j < b; ++j) {
                const auto alpha = need_a ? contract(va, j) : DensityDerivatives{};
                DensityDerivatives beta;
                if (beta_is_alpha)
                    beta = alpha;
                else if (need_b)
                    beta = contract(vb, j);
                const auto slot = range.start + static_cast<std::size_t>(j);
                store(result, DensityComponent::Alpha, slot, alpha);
                store(result, DensityComponent::Beta, slot, beta);
                if (result.contains(DensityComponent::Total))
                    store(result, DensityComponent::Total, slot, combine(alpha, beta, 1.0, 1.0));
                if (result.contains(DensityComponent::Difference))
                    store(result, DensityComponent::Difference, slot,
                          combine(alpha, beta, -1.0, 1.0));
            }
        });
    return result;
}

Field Molecule::density_pointwise(const Grid& grid, DensityComponent spin) const {
    const auto count = detail::validated_point_count(grid);
    std::vector<double> values(count);
    Eigen::MatrixXd storage;
    const Eigen::MatrixXd& p = detail::select_density_matrix(total_, spin_, spin, storage);
    for (std::size_t q = 0; q < count; ++q) {
        const Eigen::VectorXd ao = evaluator_.evaluate(detail::grid_point(grid, q));
        values[q] = ao.dot(p * ao);
    }
    return Field(grid, atoms_, std::move(values));
}

Field Molecule::density(const Grid& grid, DensityComponent spin) const {
    auto fields = densities_blocked(grid, {spin});
    return std::move(fields.at(spin));
}

SpinDensities Molecule::spin_densities(const Grid& grid) const {
    auto fields = densities_blocked(grid, {DensityComponent::Total, DensityComponent::Alpha,
                                           DensityComponent::Beta, DensityComponent::Difference});
    return {std::move(fields.at(DensityComponent::Total)),
            std::move(fields.at(DensityComponent::Alpha)),
            std::move(fields.at(DensityComponent::Beta)),
            std::move(fields.at(DensityComponent::Difference))};
}

std::map<DensityComponent, Field>
Molecule::densities(const Grid& grid, std::initializer_list<DensityComponent> components) const {
    return densities_blocked(grid, components);
}

std::map<DensityComponent, Field>
Molecule::densities_blocked(const Grid& grid, std::initializer_list<DensityComponent> components,
                            std::size_t block_size, double screening_threshold,
                            const ProgressCallback& progress) const {
    const auto count = detail::validated_point_count(grid);
    if (block_size == 0 || block_size > max_block_size)
        throw std::invalid_argument("block_size must be in [1, " + std::to_string(max_block_size) +
                                    "]");
    detail::validate_screening_threshold(screening_threshold);
    std::map<DensityComponent, std::vector<double>> buffers;
    for (DensityComponent component : components)
        buffers.try_emplace(component, count);
    if (buffers.empty()) return {};
    // Each component is a quadratic form phi^T P phi (phi = AO values at a
    // point), and each P costs one nAO x nAO x points matrix product
    // (nAO = number of AOs), which dominates the run time. Total and
    // spin difference alone use the FCHK matrices directly, one product each.
    // If alpha or beta is requested, contract the alpha and beta matrices and
    // derive total = alpha + beta and difference = alpha - beta from them, so
    // that no request needs more than two products.
    const bool by_spin =
        buffers.contains(DensityComponent::Alpha) || buffers.contains(DensityComponent::Beta);
    std::vector<DensityComponent> contracted;
    if (by_spin) {
        const bool derived = buffers.contains(DensityComponent::Total) ||
                             buffers.contains(DensityComponent::Difference);
        if (derived || buffers.contains(DensityComponent::Alpha))
            contracted.push_back(DensityComponent::Alpha);
        if (derived || buffers.contains(DensityComponent::Beta))
            contracted.push_back(DensityComponent::Beta);
    } else {
        for (const auto& entry : buffers)
            contracted.push_back(entry.first);
    }
    std::vector<Eigen::MatrixXd> storage(contracted.size());
    std::vector<const Eigen::MatrixXd*> matrices;
    for (std::size_t i = 0; i < contracted.size(); ++i)
        matrices.push_back(
            &detail::select_density_matrix(total_, spin_, contracted[i], storage[i]));
    const Eigen::Index m = static_cast<Eigen::Index>(basis_size());
    // Each block owns its AO panels and writes to disjoint output indices of
    // the preallocated buffers (the map itself is only read in the loop).
    // Keep the same arithmetic and point ordering in serial and OpenMP builds.
    detail::Progress report(progress, count);
    detail::for_each_point_block(
        count, block_size, "too many grid blocks for parallel evaluation",
        [&](detail::PointBlock range) {
            const auto [start, size] = range;
            const auto b = static_cast<Eigen::Index>(size);
            const auto points = detail::grid_points(grid, start, size);
            Eigen::MatrixXd phi(m, b);
            evaluator_.evaluate_batch(points, phi);
            // Contract only the AOs that are significant somewhere in the
            // block; if all of them are, use the full matrices as they are.
            const auto kept = detail::significant_aos(phi, screening_threshold);
            const bool screened = static_cast<Eigen::Index>(kept.size()) < m;
            if (screened) phi = detail::select_rows(phi, kept);
            std::map<DensityComponent, Eigen::VectorXd> rho;
            for (std::size_t i = 0; i < contracted.size(); ++i) {
                const Eigen::MatrixXd y =
                    screened ? Eigen::MatrixXd(detail::select_submatrix(*matrices[i], kept) * phi)
                             : Eigen::MatrixXd(*matrices[i] * phi);
                rho[contracted[i]] = (phi.array() * y.array()).colwise().sum().transpose();
            }
            for (auto& [spin, values] : buffers) {
                Eigen::Map<Eigen::VectorXd> out(values.data() + start, b);
                if (rho.contains(spin))
                    out = rho.at(spin);
                else if (spin == DensityComponent::Total)
                    out = rho.at(DensityComponent::Alpha) + rho.at(DensityComponent::Beta);
                else
                    out = rho.at(DensityComponent::Alpha) - rho.at(DensityComponent::Beta);
            }
            report.advance(static_cast<std::size_t>(b));
        });
    std::map<DensityComponent, Field> result;
    for (auto& [spin, values] : buffers)
        result.emplace(spin, Field(grid, atoms_, std::move(values)));
    return result;
}
} // namespace chem::kumo
