#include "chem/kumo/density.hpp"
#include "derivative_layout.hpp"
#include "grid_geometry.hpp"
#include "parallel_blocks.hpp"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace chem::kumo {

// The MO amplitude is a linear contraction of FCHK coefficients with AO
// values. Keep the signed result: squaring would destroy nodal information.

namespace {
void check_channel(Spin s) {
    if (s != Spin::Alpha && s != Spin::Beta)
        throw std::invalid_argument("MO spin channel must be alpha or beta");
}
} // namespace

// Restricted calculations have one set of MOs for both channels.
const Eigen::MatrixXd& Molecule::mo_coefficients(Spin channel) const {
    check_channel(channel);
    return (channel == Spin::Beta && unrestricted_) ? beta_mos_ : alpha_mos_;
}

std::size_t Molecule::orbital_count(Spin channel) const {
    return static_cast<std::size_t>(mo_coefficients(channel).cols());
}

int Molecule::occupied_count(Spin channel) const {
    check_channel(channel);
    return channel == Spin::Alpha ? nalpha_ : nbeta_;
}

const std::vector<double>& Molecule::orbital_energies(Spin channel) const {
    check_channel(channel);
    return (channel == Spin::Beta && unrestricted_) ? beta_energies_ : alpha_energies_;
}

OrbitalDerivatives Molecule::orbital_derivatives(const Eigen::Vector3d& point, std::size_t index,
                                                 Spin channel) const {
    const auto& coeff = mo_coefficients(channel);
    if (!point.allFinite()) throw std::invalid_argument("MO evaluation point must be finite");
    if (index == 0 || index > static_cast<std::size_t>(coeff.cols()))
        throw std::out_of_range("MO index outside available coefficient range");
    const auto ao = evaluator_.evaluate_derivatives(point);
    const auto c = coeff.col(static_cast<Eigen::Index>(index - 1));
    OrbitalDerivatives result;
    result.value = c.dot(ao.values);
    result.gradient = ao.gradient.transpose() * c;
    const Eigen::Matrix<double, 6, 1> h = ao.hessian.transpose() * c;
    result.hessian(0, 0) = h[0];
    result.hessian(1, 1) = h[1];
    result.hessian(2, 2) = h[2];
    result.hessian(0, 1) = result.hessian(1, 0) = h[3];
    result.hessian(0, 2) = result.hessian(2, 0) = h[4];
    result.hessian(1, 2) = result.hessian(2, 1) = h[5];
    return result;
}

std::vector<OrbitalDerivatives> Molecule::orbital_derivatives_blocked(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points, std::size_t index,
    Spin channel, std::size_t block_size) const {
    const auto& coeff = mo_coefficients(channel);
    if (index == 0 || index > static_cast<std::size_t>(coeff.cols()))
        throw std::out_of_range("MO index outside available coefficient range");
    if (block_size == 0)
        throw std::invalid_argument("orbital derivative block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("MO evaluation points must be finite");
    const std::size_t count = static_cast<std::size_t>(points.cols());
    std::vector<OrbitalDerivatives> result(count);
    const auto c = coeff.col(static_cast<Eigen::Index>(index - 1));
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto capacity = static_cast<Eigen::Index>(std::min(block_size, count));
    // Each worker reuses its AO panels: value, x, y, z, xx, yy, zz, xy, xz,
    // yz (one column per point each).
    detail::for_each_point_block_with_workspace(
        count, block_size, "too many orbital derivative blocks for parallel evaluation",
        [&] { return Eigen::MatrixXd(n, 10 * capacity); },
        [&](detail::PointBlock range, Eigen::MatrixXd& ao) {
            const auto b = static_cast<Eigen::Index>(range.size);
            evaluator_.evaluate_derivatives_batch(
                points.middleCols(static_cast<Eigen::Index>(range.start), b), 2,
                ao.leftCols(10 * b));
            const auto panel = [&](int component, Eigen::Index j) {
                return c.dot(ao.col(component * b + j));
            };
            for (Eigen::Index j = 0; j < b; ++j) {
                auto& out = result[range.start + static_cast<std::size_t>(j)];
                out.value = panel(0, j);
                for (int a = 0; a < 3; ++a) {
                    out.gradient[a] = panel(1 + a, j);
                    for (int e = a; e < 3; ++e)
                        out.hessian(a, e) = out.hessian(e, a) =
                            panel(4 + detail::ao_hessian_column[a][e], j);
                }
            }
        });
    return result;
}

Field Molecule::orbital(const Grid& grid, std::size_t index, Spin channel) const {
    auto fields = orbitals(grid, {index}, channel);
    return std::move(fields.front());
}

std::vector<Field> Molecule::orbitals(const Grid& grid, const std::vector<std::size_t>& indices,
                                      Spin channel, const ProgressCallback& progress) const {
    const auto& coeff = mo_coefficients(channel);
    if (coeff.cols() == 0) throw std::runtime_error("FCHK has no MO coefficients");
    for (auto i : indices)
        if (i == 0 || i > static_cast<std::size_t>(coeff.cols()))
            throw std::out_of_range("MO index outside available coefficient range");
    if (indices.empty()) return {};
    const auto count = detail::validated_point_count(grid);
    std::vector<std::vector<double>> buffers(indices.size(), std::vector<double>(count));
    Eigen::MatrixXd selected(coeff.rows(), static_cast<Eigen::Index>(indices.size()));
    for (std::size_t i = 0; i < indices.size(); ++i)
        selected.col(static_cast<Eigen::Index>(i)) =
            coeff.col(static_cast<Eigen::Index>(indices[i] - 1));
    constexpr std::size_t block_size = 128;
    detail::Progress report(progress, count);
    detail::for_each_point_block(
        count, block_size, "too many orbital blocks for parallel evaluation",
        [&](detail::PointBlock range) {
            const auto [start, size] = range;
            const auto b = static_cast<Eigen::Index>(size);
            const auto points = detail::grid_points(grid, start, size);
            Eigen::MatrixXd phi(static_cast<Eigen::Index>(basis_size()), b);
            evaluator_.evaluate_batch(points, phi);
            const Eigen::MatrixXd psi = selected.transpose() * phi;
            for (std::size_t mo = 0; mo < indices.size(); ++mo)
                for (Eigen::Index t = 0; t < b; ++t)
                    buffers[mo][start + static_cast<std::size_t>(t)] =
                        psi(static_cast<Eigen::Index>(mo), t);
            report.advance(size);
        });
    std::vector<Field> result;
    result.reserve(indices.size());
    for (auto& buffer : buffers)
        result.emplace_back(grid, atoms_, std::move(buffer));
    return result;
}
} // namespace chem::kumo
