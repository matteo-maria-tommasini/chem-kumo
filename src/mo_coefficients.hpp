#pragma once

#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace chem::kumo::detail {

/// Read an MO coefficient record ("Alpha MO coefficients" or
/// "Beta MO coefficients") as an nbasis x nmo matrix, one MO per column.
///
/// Gaussian FCHK stores the coefficients MO by MO, with the AO index running
/// fastest. The record must exist; callers decide whether a missing record
/// is an error. Throws std::runtime_error if the record size is not a whole
/// number of MOs, if it has more MOs than basis functions, or if any
/// coefficient is not finite.
inline Eigen::MatrixXd read_mo_coefficients(const FormattedCheckpoint& fchk,
                                            const std::string& label, std::size_t nbasis) {
    const auto& flat = fchk.reals(label);
    if (nbasis == 0 || flat.size() % nbasis != 0 || flat.size() / nbasis > nbasis)
        throw std::runtime_error("invalid MO coefficient dimensions: " + label);
    const std::size_t nmo = flat.size() / nbasis;
    Eigen::MatrixXd coefficients(static_cast<Eigen::Index>(nbasis), static_cast<Eigen::Index>(nmo));
    for (std::size_t mo = 0; mo < nmo; ++mo)
        for (std::size_t ao = 0; ao < nbasis; ++ao)
            coefficients(static_cast<Eigen::Index>(ao), static_cast<Eigen::Index>(mo)) =
                flat[mo * nbasis + ao];
    if (!coefficients.allFinite()) throw std::runtime_error("nonfinite MO coefficients: " + label);
    return coefficients;
}

} // namespace chem::kumo::detail
