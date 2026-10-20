#pragma once

#include "chem/kumo/local_shell_integrals.hpp"

#include <Eigen/Core>
#include <cstddef>
#include <vector>

namespace chem::kumo::detail {

// Accumulate one canonical shell-pair quartet into the upper AO triangle.
// Keep the established shell/AO multiplicity rules unchanged.
// The caller is responsible for canonical pair traversal and density validation.
inline void accumulate_local_density_quartet(
    Eigen::MatrixXd& upper,
    const Eigen::Ref<const Eigen::MatrixXd>& density,
    const std::vector<Eigen::Index>& offsets,
    std::size_t a, std::size_t b, std::size_t c, std::size_t d,
    bool same_shell_pair,
    const LocalShellQuartet& q) {
    const double mult_ab = a == b ? 1.0 : 2.0;
    const double mult_cd = c == d ? 1.0 : 2.0;
    for (std::size_t i = 0; i < q.dimensions[0]; ++i)
        for (std::size_t j = 0; j < q.dimensions[1]; ++j) {
            const Eigen::Index mu = offsets[a] + static_cast<Eigen::Index>(i);
            const Eigen::Index nu = offsets[b] + static_cast<Eigen::Index>(j);
            if (mu > nu) continue;
            for (std::size_t k = 0; k < q.dimensions[2]; ++k)
                for (std::size_t l = 0; l < q.dimensions[3]; ++l) {
                    const Eigen::Index lambda = offsets[c] + static_cast<Eigen::Index>(k);
                    const Eigen::Index sigma = offsets[d] + static_cast<Eigen::Index>(l);
                    const double integral = q.at(i, j, k, l);
                    upper(mu, nu) += mult_cd * density(lambda, sigma) * integral;
                    if (!same_shell_pair && lambda <= sigma)
                        upper(lambda, sigma) += ((a == b && mu != nu) ? 2.0 : mult_ab)
                                                * density(mu, nu) * integral;
                }
        }
}

} // namespace chem::kumo::detail
