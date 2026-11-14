#pragma once

// Test-only access to the private evaluation kernels and data of Molecule.
// Since chem-kumo 2.0 these are not public: applications use Evaluator. The
// tests still check each kernel on its own, and use the simple pointwise
// implementations as independent references.

#include "chem/kumo/density.hpp"

#include <Eigen/Core>
#include <initializer_list>
#include <map>
#include <utility>
#include <vector>

namespace chem::kumo::detail {
struct MoleculeTestAccess {
    // Density kernels from the density matrices (AO_P).

    static Field density(const Molecule& m, const Grid& grid,
                         DensityComponent component = DensityComponent::Total) {
        return m.densities_blocked(grid, {component}).at(component);
    }

    static std::map<DensityComponent, Field> densities_blocked(
        const Molecule& m, const Grid& grid, std::initializer_list<DensityComponent> components,
        std::size_t block_size = 128, double screening_threshold = default_screening_threshold,
        const ProgressCallback& progress = {}) {
        return m.densities_blocked(grid, components, block_size, screening_threshold, progress);
    }

    /// The four density components on one grid, as the removed
    /// Molecule::spin_densities() returned them.
    struct SpinDensities {
        Field total, alpha, beta, spin;
    };

    static SpinDensities spin_densities(const Molecule& m, const Grid& grid) {
        auto f = m.densities_blocked(grid, {DensityComponent::Total, DensityComponent::Alpha,
                                            DensityComponent::Beta, DensityComponent::Difference});
        return {std::move(f.at(DensityComponent::Total)), std::move(f.at(DensityComponent::Alpha)),
                std::move(f.at(DensityComponent::Beta)),
                std::move(f.at(DensityComponent::Difference))};
    }

    static Field density_pointwise(const Molecule& m, const Grid& grid,
                                   DensityComponent component = DensityComponent::Total) {
        return m.density_pointwise(grid, component);
    }

    static DensityDerivatives
    density_derivatives(const Molecule& m, const Eigen::Vector3d& point,
                        DensityComponent component = DensityComponent::Total) {
        return m.density_derivatives(point, component);
    }

    static std::vector<DensityDerivatives> density_derivatives_blocked(
        const Molecule& m, const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        DensityComponent component = DensityComponent::Total, std::size_t block_size = 128,
        double screening_threshold = default_screening_threshold) {
        return m.density_derivatives_blocked(points, component, block_size, screening_threshold);
    }

    static std::map<DensityComponent, std::vector<DensityDerivatives>> density_derivatives_shared(
        const Molecule& m, const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        const std::vector<DensityComponent>& components, std::size_t block_size, int ao_order = 2,
        double screening_threshold = default_screening_threshold) {
        return m.density_derivatives_shared(points, components, block_size, ao_order,
                                            screening_threshold);
    }

    // Orbital kernels.

    static Field orbital(const Molecule& m, const Grid& grid, std::size_t index,
                         SpinChannel channel = SpinChannel::Alpha) {
        return m.orbitals(grid, {index}, channel).front();
    }

    static std::vector<Field> orbitals(const Molecule& m, const Grid& grid,
                                       const std::vector<std::size_t>& indices,
                                       SpinChannel channel = SpinChannel::Alpha,
                                       const ProgressCallback& progress = {}) {
        return m.orbitals(grid, indices, channel, progress);
    }

    static OrbitalDerivatives orbital_derivatives(const Molecule& m, const Eigen::Vector3d& point,
                                                  std::size_t index,
                                                  SpinChannel channel = SpinChannel::Alpha) {
        return m.orbital_derivatives(point, index, channel);
    }

    // AO_MO kernel and the occupied-orbital density matrices.

    static auto density_derivatives_mo(
        const Molecule& m, const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        const std::vector<DensityComponent>& components, std::size_t block, int order) {
        return m.density_derivatives_mo(points, components, block, order);
    }

    static auto reconstructed_matrices(const Molecule& m) {
        const auto& cb = m.unrestricted_ ? m.beta_mos_ : m.alpha_mos_;
        const Eigen::MatrixXd pa =
            m.alpha_mos_.leftCols(m.nalpha_) * m.alpha_mos_.leftCols(m.nalpha_).transpose();
        const Eigen::MatrixXd pb = cb.leftCols(m.nbeta_) * cb.leftCols(m.nbeta_).transpose();
        std::map<DensityComponent, Eigen::MatrixXd> result;
        result.emplace(DensityComponent::Alpha, pa);
        result.emplace(DensityComponent::Beta, pb);
        result.emplace(DensityComponent::Total, pa + pb);
        result.emplace(DensityComponent::Difference, pa - pb);
        return result;
    }

    static DensityDerivatives ao_p_reference(const Molecule& m, const Eigen::Vector3d& point,
                                             const Eigen::MatrixXd& p) {
        const auto ao = m.evaluator_.evaluate_derivatives(point, 2);
        const auto y = p * ao.values;
        const auto pg = p * ao.gradient;
        constexpr int index[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
        DensityDerivatives d;
        d.value = ao.values.dot(y);
        for (int a = 0; a < 3; ++a) {
            d.gradient[a] = 2.0 * ao.gradient.col(a).dot(y);
            d.tau += 0.5 * ao.gradient.col(a).dot(pg.col(a));
            for (int b = a; b < 3; ++b)
                d.hessian(a, b) = d.hessian(b, a) =
                    2.0 * (ao.hessian.col(index[a][b]).dot(y) + ao.gradient.col(a).dot(pg.col(b)));
        }
        d.laplacian = d.hessian.trace();
        return d;
    }
};
} // namespace chem::kumo::detail

namespace kumo_test {
/// Short name for the test access to the Molecule kernels.
using Kernels = chem::kumo::detail::MoleculeTestAccess;
} // namespace kumo_test
