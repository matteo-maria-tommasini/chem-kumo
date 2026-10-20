#pragma once

#include "chem/kumo/local_shell_integrals.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

/// Exact L2 norms of AO products, B(mu,nu) = sqrt(int (chi_mu chi_nu)^2 dr).
/// Diagnostic O(number of shell pairs) preprocessing; no four-index tensor.
[[nodiscard]] inline Eigen::MatrixXd local_ao_product_norms(const BasisSet& basis) {
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(offsets.back(), offsets.back());
    for (std::size_t a = 0; a < shells.size(); ++a)
        for (std::size_t b = a; b < shells.size(); ++b) {
            const auto q = local_shell_quartet({&shells[a], &shells[b], &shells[a], &shells[b]});
            for (std::size_t i = 0; i < shells[a].size(); ++i)
                for (std::size_t j = 0; j < shells[b].size(); ++j) {
                    const double squared = q.at(i, j, i, j);
                    if (!std::isfinite(squared) || squared < -1e-12 * std::max(1.0, std::abs(squared)))
                        throw std::runtime_error("invalid AO product norm squared");
                    const double value = std::sqrt(std::max(0.0, squared));
                    const Eigen::Index mu = offsets[a] + static_cast<Eigen::Index>(i);
                    const Eigen::Index nu = offsets[b] + static_cast<Eigen::Index>(j);
                    norms(mu, nu) = norms(nu, mu) = value;
                }
        }
    return norms;
}

/// Upper bound for one omitted ordered AO-density contribution.
[[nodiscard]] inline double local_ao_contribution_bound(
    const Eigen::Ref<const Eigen::MatrixXd>& norms,
    const Eigen::Ref<const Eigen::MatrixXd>& density,
    Eigen::Index mu, Eigen::Index nu, Eigen::Index lambda, Eigen::Index sigma) {
    if (norms.rows() != norms.cols() || density.rows() != norms.rows() ||
        density.cols() != norms.cols())
        throw std::invalid_argument("AO norm/density shape mismatch");
    return std::abs(density(lambda, sigma)) * norms(mu, nu) * norms(lambda, sigma);
}

} // namespace chem::kumo
