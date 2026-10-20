#pragma once

#include "chem/kumo/local_density_quartet_contraction.hpp"

#include "chem/kumo/local_shell_pair_integrals.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chem::kumo {

/// Unscreened symmetric density contraction with reusable primitive-pair data.
/// All canonical shell pairs are prepared once; no four-index tensor is stored.
/// This reference implementation is intended for small and medium bases.
[[nodiscard]] inline Eigen::MatrixXd local_density_overlap_symmetric_paired(
    const BasisSet& basis, const Eigen::Ref<const Eigen::MatrixXd>& density) {
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao = offsets.back();
    if (density.rows() != nao || density.cols() != nao)
        throw std::invalid_argument("local_density_overlap_symmetric_paired: AO density shape mismatch");
    if (!density.allFinite())
        throw std::invalid_argument("local_density_overlap_symmetric_paired: nonfinite AO density");
    if (nao == 0) return Eigen::MatrixXd(0, 0);
    const double asym = (density - density.transpose()).cwiseAbs().maxCoeff();
    const double scale = density.cwiseAbs().maxCoeff();
    if (asym > 1e-10 * std::max(1.0, scale))
        throw std::invalid_argument("local_density_overlap_symmetric_paired: asymmetric AO density");

    using Pair = std::pair<std::size_t, std::size_t>;
    std::vector<Pair> pairs;
    for (std::size_t a = 0; a < shells.size(); ++a)
        for (std::size_t b = a; b < shells.size(); ++b)
            pairs.emplace_back(a, b);
    std::vector<LocalPreparedShellPair> prepared;
    prepared.reserve(pairs.size());
    for (const auto& [a, b] : pairs)
        prepared.push_back(prepare_local_shell_pair(shells[a], shells[b]));

    Eigen::MatrixXd upper = Eigen::MatrixXd::Zero(nao, nao);
    for (std::size_t ip = 0; ip < pairs.size(); ++ip) {
        const auto [a, b] = pairs[ip];
        for (std::size_t iq = ip; iq < pairs.size(); ++iq) {
            const auto [c, d] = pairs[iq];
            const auto q = local_shell_quartet_paired(prepared[ip], prepared[iq]);
            detail::accumulate_local_density_quartet(upper, density, offsets, a, b, c, d, ip == iq, q);
        }
    }
    return Eigen::MatrixXd(upper.selfadjointView<Eigen::Upper>());
}

} // namespace chem::kumo
