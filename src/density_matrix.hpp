#pragma once

// AO density matrices of the four spin components, from the two canonical
// FCHK matrices (total and spin difference).

#include "chem/kumo/density.hpp"

#include <Eigen/Core>
#include <stdexcept>

namespace chem::kumo::detail {
// AO density matrix of one spin component. Total and spin difference are the
// canonical FCHK matrices and are returned by reference; alpha and beta are
// built in `storage` as (P_total +/- P_spin) / 2, never from MO occupations.
inline const Eigen::MatrixXd& select_density_matrix(const Eigen::MatrixXd& total,
                                                    const Eigen::MatrixXd& spin, Spin component,
                                                    Eigen::MatrixXd& storage) {
    switch (component) {
    case Spin::Total: return total;
    case Spin::Difference: return spin;
    case Spin::Alpha: storage = 0.5 * (total + spin); return storage;
    case Spin::Beta: storage = 0.5 * (total - spin); return storage;
    }
    throw std::invalid_argument("invalid density spin component");
}
} // namespace chem::kumo::detail
