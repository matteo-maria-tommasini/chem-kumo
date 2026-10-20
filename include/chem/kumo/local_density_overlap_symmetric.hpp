#pragma once

#include "chem/kumo/local_density_quartet_contraction.hpp"

#include "chem/kumo/local_shell_integrals.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chem::kumo {

/// Unscreened, symmetry-reduced contraction of local four-AO integrals.
/// Uses unique unordered shell pairs and pair-of-pairs; no global ERI tensor.
/// P is the symmetric contravariant AO density (rho = chi^T P chi).
[[nodiscard]] inline Eigen::MatrixXd local_density_overlap_symmetric(
    const BasisSet& basis, const Eigen::Ref<const Eigen::MatrixXd>& density) {
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao = offsets.back();
    if (density.rows() != nao || density.cols() != nao)
        throw std::invalid_argument("local_density_overlap_symmetric: AO density shape mismatch");
    if (!density.allFinite())
        throw std::invalid_argument("local_density_overlap_symmetric: nonfinite AO density");
    if (nao == 0) return Eigen::MatrixXd(0, 0);
    const double asym = (density - density.transpose()).cwiseAbs().maxCoeff();
    const double scale = density.cwiseAbs().maxCoeff();
    if (asym > 1e-10 * std::max(1.0, scale))
        throw std::invalid_argument("local_density_overlap_symmetric: asymmetric AO density");

    using Pair = std::pair<std::size_t, std::size_t>;
    std::vector<Pair> pairs;
    for (std::size_t a = 0; a < shells.size(); ++a)
        for (std::size_t b = a; b < shells.size(); ++b)
            pairs.emplace_back(a, b);

    // Only upper-triangular shell blocks are accumulated. Mirror at the end.
    Eigen::MatrixXd upper = Eigen::MatrixXd::Zero(nao, nao);
    for (std::size_t ip = 0; ip < pairs.size(); ++ip) {
        const auto [a, b] = pairs[ip];
        for (std::size_t iq = ip; iq < pairs.size(); ++iq) {
            const auto [c, d] = pairs[iq];
            const auto q = local_shell_quartet({&shells[a], &shells[b], &shells[c], &shells[d]});
            detail::accumulate_local_density_quartet(upper, density, offsets, a, b, c, d, ip == iq, q);
        }
    }
    Eigen::MatrixXd result = upper.selfadjointView<Eigen::Upper>();
    return result;
}

} // namespace chem::kumo
