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

// A molecular-orbital (MO) amplitude is the sum of atomic-orbital (AO) values
// weighted by the FCHK MO coefficients. The signed value is returned:
// squaring it would lose the nodal structure of the orbital.

namespace {
void check_channel(SpinChannel s) {
    if (s != SpinChannel::Alpha && s != SpinChannel::Beta)
        throw std::invalid_argument("MO spin channel must be alpha or beta");
}
} // namespace

// Restricted calculations have one set of MOs for both channels.
const Eigen::MatrixXd& Molecule::mo_coefficients(SpinChannel channel) const {
    check_channel(channel);
    return (channel == SpinChannel::Beta && unrestricted_) ? beta_mos_ : alpha_mos_;
}

std::size_t Molecule::orbital_count(SpinChannel channel) const {
    return static_cast<std::size_t>(mo_coefficients(channel).cols());
}

std::size_t Molecule::occupied_count(SpinChannel channel) const {
    check_channel(channel);
    return static_cast<std::size_t>(channel == SpinChannel::Alpha ? nalpha_ : nbeta_);
}

const std::vector<double>& Molecule::orbital_energies(SpinChannel channel) const {
    check_channel(channel);
    return (channel == SpinChannel::Beta && unrestricted_) ? beta_energies_ : alpha_energies_;
}

detail::OrbitalDerivatives Molecule::orbital_derivatives(const Eigen::Vector3d& point,
                                                         std::size_t index,
                                                         SpinChannel channel) const {
    const auto& coeff = mo_coefficients(channel);
    if (!point.allFinite()) throw std::invalid_argument("MO evaluation point must be finite");
    if (index == 0 || index > static_cast<std::size_t>(coeff.cols()))
        throw std::out_of_range("MO index outside available coefficient range");
    const auto ao = evaluator_.evaluate_derivatives(point);
    const auto c = coeff.col(static_cast<Eigen::Index>(index - 1));
    detail::OrbitalDerivatives result;
    result.value = c.dot(ao.values);
    result.gradient = ao.gradient.transpose() * c;
    const Eigen::Matrix<double, 6, 1> h = ao.hessian.transpose() * c;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.hessian(a, b) = h[detail::ao_hessian_column[a][b]];
    return result;
}

std::vector<Field> Molecule::orbitals(const Grid& grid, const std::vector<std::size_t>& indices,
                                      SpinChannel channel, const ProgressCallback& progress) const {
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
