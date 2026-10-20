#pragma once

#include "chem/kumo/basis_evaluator.hpp"
#include "chem/kumo/local_gaussian_integrals.hpp"

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

/// An AO component in Gaussian FCHK order, expanded in normalized Cartesian primitives.
using LocalContractedComponent = std::vector<LocalGaussianPrimitive>;

/// One shell (or one part of an SP shell), in Gaussian FCHK component order.
using LocalShellComponents = std::vector<LocalContractedComponent>;

/// Convert an FCHK basis into contracted AO components. SP shells are split s then p.
/// The Cartesian-to-pure convention is shared with BasisEvaluator.
[[nodiscard]] inline std::vector<LocalShellComponents> local_shell_components(const BasisSet& basis) {
    std::vector<LocalShellComponents> result;
    const auto append = [&](const Shell& shell, int l, const std::vector<double>& coeff,
                            ShellRepresentation rep) {
        const auto powers = gaussian_cartesian_order(l);
        const Eigen::Index nc = static_cast<Eigen::Index>(powers.size());
        Eigen::MatrixXd transform = Eigen::MatrixXd::Identity(nc, nc);
        if (rep == ShellRepresentation::Pure && l >= 2)
            transform = gaussian_cartesian_to_pure(l);
        LocalShellComponents components(static_cast<std::size_t>(transform.rows()));
        if (shell.exponents.size() != coeff.size())
            throw std::invalid_argument("shell exponent/coefficient count mismatch");
        for (Eigen::Index out = 0; out < transform.rows(); ++out)
            for (Eigen::Index cart = 0; cart < nc; ++cart) {
                const double factor = transform(out, cart);
                if (factor == 0.0) continue;
                const auto& power = powers[static_cast<std::size_t>(cart)];
                const int L = power.x + power.y + power.z;
                const auto odd = [](int n) {
                    double value = 1.0;
                    for (int i = n; i > 0; i -= 2) value *= i;
                    return value;
                };
                for (std::size_t p = 0; p < shell.exponents.size(); ++p) {
                    const double alpha = shell.exponents[p];
                    if (!(alpha > 0.0)) throw std::invalid_argument("nonpositive exponent");
                    const double norm = std::pow(2.0 * alpha / std::numbers::pi, 0.75) *
                        std::sqrt(std::pow(4.0 * alpha, L) /
                            (odd(2 * power.x - 1) * odd(2 * power.y - 1) *
                             odd(2 * power.z - 1)));
                    components[static_cast<std::size_t>(out)].push_back({
                        {shell.origin.x(), shell.origin.y(), shell.origin.z()},
                        {power.x, power.y, power.z}, alpha, factor * coeff[p] * norm});
                }
            }
        result.push_back(std::move(components));
    };
    for (const auto& shell : basis.shells()) {
        if (shell.is_sp) {
            append(shell, 0, shell.coefficients, ShellRepresentation::Cartesian);
            append(shell, 1, shell.p_coefficients, ShellRepresentation::Cartesian);
        } else {
            append(shell, shell.angular_momentum, shell.coefficients, shell.representation);
        }
    }
    return result;
}

/// Flat tensor of one contracted shell quartet, dimensions in FCHK AO order.
/// Index: ((i*n1+j)*n2+k)*n3+l. Only the requested quartet is materialized.
struct LocalShellQuartet {
    std::array<std::size_t, 4> dimensions{};
    std::vector<double> values;
    [[nodiscard]] double at(std::size_t i, std::size_t j,
                            std::size_t k, std::size_t l) const {
        if (i >= dimensions[0] || j >= dimensions[1] ||
            k >= dimensions[2] || l >= dimensions[3])
            throw std::out_of_range("shell quartet component index");
        return values[((i * dimensions[1] + j) * dimensions[2] + k) * dimensions[3] + l];
    }
};

[[nodiscard]] inline LocalShellQuartet local_shell_quartet(
    const std::array<const LocalShellComponents*, 4>& shells) {
    LocalShellQuartet result;
    std::size_t size = 1;
    for (std::size_t d = 0; d < 4; ++d) {
        if (!shells[d]) throw std::invalid_argument("null shell");
        result.dimensions[d] = shells[d]->size();
        size *= result.dimensions[d];
    }
    result.values.reserve(size);
    for (const auto& a : *shells[0])
        for (const auto& b : *shells[1])
            for (const auto& c : *shells[2])
                for (const auto& d : *shells[3]) {
                    double value = 0.0;
                    for (const auto& pa : a)
                        for (const auto& pb : b)
                            for (const auto& pc : c)
                                for (const auto& pd : d)
                                    value += local_four_gaussian_integral({pa, pb, pc, pd});
                    result.values.push_back(value);
                }
    return result;
}

} // namespace chem::kumo
