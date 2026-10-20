#pragma once

#include "chem/kumo/local_shell_integrals.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chem::kumo {

/// Exact local four-AO contraction with an AO density matrix.
/// The matrix P is the contravariant AO density (rho = chi^T P chi).
/// No global four-index tensor is formed. This unscreened reference
/// implementation is O(N_shell^4) and is intended for small systems.
[[nodiscard]] inline Eigen::MatrixXd local_density_overlap(
    const BasisSet& basis, const Eigen::Ref<const Eigen::MatrixXd>& density) {
    const auto shells = local_shell_components(basis);
    std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
    for (std::size_t a = 0; a < shells.size(); ++a)
        offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
    const Eigen::Index nao = offsets.back();
    if (density.rows() != nao || density.cols() != nao)
        throw std::invalid_argument("local_density_overlap: AO density shape mismatch");
    if (!density.allFinite())
        throw std::invalid_argument("local_density_overlap: nonfinite AO density");
    if (nao == 0) return Eigen::MatrixXd(0, 0);
    const double asym = (density - density.transpose()).cwiseAbs().maxCoeff();
    const double scale = density.cwiseAbs().maxCoeff();
    if (asym > 1e-10 * std::max(1.0, scale))
        throw std::invalid_argument("local_density_overlap: asymmetric AO density");

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nao, nao);
    // Only the external shell-pair symmetry is used. Internal pairs
    // (c,d) are traversed in full, so off-diagonal P contributions
    // have their correct multiplicity without special-case factors.
    for (std::size_t a = 0; a < shells.size(); ++a)
        for (std::size_t b = a; b < shells.size(); ++b) {
            Eigen::MatrixXd block = Eigen::MatrixXd::Zero(
                static_cast<Eigen::Index>(shells[a].size()),
                static_cast<Eigen::Index>(shells[b].size()));
            for (std::size_t c = 0; c < shells.size(); ++c)
                for (std::size_t d = 0; d < shells.size(); ++d) {
                    const auto q = local_shell_quartet({&shells[a], &shells[b], &shells[c], &shells[d]});
                    for (std::size_t i = 0; i < q.dimensions[0]; ++i)
                        for (std::size_t j = 0; j < q.dimensions[1]; ++j) {
                            double contribution = 0.0;
                            for (std::size_t k = 0; k < q.dimensions[2]; ++k)
                                for (std::size_t l = 0; l < q.dimensions[3]; ++l)
                                    contribution += density(offsets[c] + static_cast<Eigen::Index>(k),
                                                            offsets[d] + static_cast<Eigen::Index>(l)) * q.at(i,j,k,l);
                            block(static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j)) += contribution;
                        }
                }
            for (Eigen::Index i = 0; i < block.rows(); ++i)
                for (Eigen::Index j = 0; j < block.cols(); ++j) {
                    const Eigen::Index mu = offsets[a] + i;
                    const Eigen::Index nu = offsets[b] + j;
                    result(mu, nu) = block(i,j);
                    result(nu, mu) = block(i,j);
                }
        }
    return result;
}

} // namespace chem::kumo
