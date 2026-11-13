#include "basis_math.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace chem::kumo::detail {

double double_factorial_odd(int n) noexcept {
    if (n <= 0) return 1.0;
    double value = 1.0;
    for (int i = n; i >= 1; i -= 2)
        value *= static_cast<double>(i);
    return value;
}

double cartesian_primitive_normalization(double alpha, const CartesianPower& p) {
    const int l = p.x + p.y + p.z;
    const double radial = std::pow(2.0 * alpha / std::numbers::pi, 0.75);
    const double angular =
        std::sqrt(std::pow(4.0 * alpha, static_cast<double>(l)) /
                  (double_factorial_odd(2 * p.x - 1) * double_factorial_odd(2 * p.y - 1) *
                   double_factorial_odd(2 * p.z - 1)));
    return radial * angular;
}

std::vector<NormalizedBlock> normalized_blocks(const BasisSet& basis) {
    std::vector<NormalizedBlock> blocks;
    blocks.reserve(basis.shells().size() + 4U);
    std::size_t offset = 0;
    const auto add = [&](const Shell& shell, int l, const std::vector<double>& coefficients,
                         ShellRepresentation representation) {
        NormalizedBlock block;
        block.origin = shell.origin;
        block.angular_momentum = l;
        block.exponents = shell.exponents;
        block.powers = gaussian_cartesian_order(l);
        const auto ncart = block.powers.size();
        const auto nprim = block.exponents.size();
        block.weights.resize(static_cast<Eigen::Index>(ncart), static_cast<Eigen::Index>(nprim));
        for (std::size_t k = 0; k < ncart; ++k)
            for (std::size_t p = 0; p < nprim; ++p)
                block.weights(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(p)) =
                    coefficients[p] *
                    cartesian_primitive_normalization(block.exponents[p], block.powers[k]);
        if (representation == ShellRepresentation::Pure && l >= 2)
            block.pure_transform = gaussian_cartesian_to_pure(l);
        block.offset = offset;
        block.count = block.pure_transform.size()
                          ? static_cast<std::size_t>(block.pure_transform.rows())
                          : ncart;
        offset += block.count;
        blocks.push_back(std::move(block));
    };
    for (const auto& shell : basis.shells()) {
        if (shell.is_sp) {
            add(shell, 0, shell.coefficients, ShellRepresentation::Cartesian);
            add(shell, 1, shell.p_coefficients, ShellRepresentation::Cartesian);
        } else {
            add(shell, shell.angular_momentum, shell.coefficients, shell.representation);
        }
    }
    if (offset != basis.function_count())
        throw std::logic_error("AO blocks do not match the basis size");
    return blocks;
}

} // namespace chem::kumo::detail
