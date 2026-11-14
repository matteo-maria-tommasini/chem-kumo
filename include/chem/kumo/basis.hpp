#pragma once

#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chem::kumo {

/// Angular form of a Gaussian shell: Cartesian (the (l + 1)(l + 2) / 2
/// monomials x^a y^b z^c with a + b + c = l) or pure (the 2l + 1 real
/// spherical harmonics).
enum class ShellRepresentation {
    Cartesian,
    Pure,
};

/// One contracted Gaussian shell as read from a Gaussian formatted
/// checkpoint file (FCHK), in the order of the file. An SP shell combines an
/// s and a p shell that share their exponents; its p contraction
/// coefficients are stored in p_coefficients.
struct Shell {
    std::size_t center = 0;   ///< zero-based atom index
    int angular_momentum = 0; ///< l (0 = s, 1 = p, ...); 1 for an SP shell
    ShellRepresentation representation = ShellRepresentation::Cartesian;
    bool is_sp = false;                               ///< SP shell (s and p, shared exponents)
    Eigen::Vector3d origin = Eigen::Vector3d::Zero(); ///< shell centre, bohr
    std::vector<double> exponents;                    ///< primitive exponents, bohr^-2
    /// Contraction coefficients of the shell (of its s part for an SP
    /// shell), multiplying normalized primitives.
    std::vector<double> coefficients;
    std::vector<double> p_coefficients; ///< p part of an SP shell, else empty

    /// Number of atomic orbitals (AOs, the basis functions) of the shell:
    /// 4 for SP, 2l + 1 for pure shells with l >= 2, (l + 1)(l + 2) / 2
    /// otherwise (s and p shells are the same in both forms).
    [[nodiscard]] std::size_t function_count() const;
};

/// Gaussian basis set rebuilt from the shell records of an FCHK, in the
/// order of the file.
class BasisSet {
public:
    /// Build the shells from the FCHK shell records. Throws
    /// std::runtime_error for missing, inconsistent or unsupported records,
    /// or if the number of AOs differs from "Number of basis functions".
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
