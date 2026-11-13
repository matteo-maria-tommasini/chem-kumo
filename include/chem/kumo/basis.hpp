#pragma once

#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chem::kumo {

/// Gaussian shell angular representation.
enum class ShellRepresentation {
    Cartesian,
    Pure,
};

/// Contracted shell metadata in the input FCHK basis ordering.
/// The optional p coefficients represent combined SP shells.
struct Shell {
    std::size_t center = 0;   ///< zero-based atom index
    int angular_momentum = 0; ///< l; 1 for the p part of an SP shell
    ShellRepresentation representation = ShellRepresentation::Cartesian;
    bool is_sp = false;                               ///< combined s and p shell
    Eigen::Vector3d origin = Eigen::Vector3d::Zero(); ///< bohr
    std::vector<double> exponents;                    ///< primitive exponents, bohr^-2
    /// Contraction coefficients of the (s part of the) shell, for normalized
    /// primitives.
    std::vector<double> coefficients;
    std::vector<double> p_coefficients; ///< p part of an SP shell, else empty

    /// Number of AOs: 4 for SP, 2l + 1 for pure shells with l >= 2,
    /// (l + 1)(l + 2) / 2 otherwise.
    [[nodiscard]] std::size_t function_count() const;
};

/// Ordered Gaussian basis reconstructed from FCHK shell records.
class BasisSet {
public:
    /// Build the shells from the FCHK records; throws std::runtime_error for
    /// missing, inconsistent or unsupported records.
    explicit BasisSet(const FormattedCheckpoint& fchk);

    /// Shells in FCHK order; their AOs are consecutive in that order.
    [[nodiscard]] const std::vector<Shell>& shells() const noexcept { return shells_; }

    /// Total number of AOs.
    [[nodiscard]] std::size_t function_count() const noexcept { return function_count_; }

private:
    std::vector<Shell> shells_;
    std::size_t function_count_ = 0;
};

} // namespace chem::kumo
