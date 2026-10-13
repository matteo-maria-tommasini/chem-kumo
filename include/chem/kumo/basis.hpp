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
    std::size_t center = 0;
    int angular_momentum = 0;
    ShellRepresentation representation = ShellRepresentation::Cartesian;
    bool is_sp = false;
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    std::vector<double> exponents;
    std::vector<double> coefficients;
    std::vector<double> p_coefficients;

    [[nodiscard]] std::size_t function_count() const;
};

/// Ordered Gaussian basis reconstructed from FCHK shell records.
class BasisSet {
public:
    explicit BasisSet(const FormattedCheckpoint& fchk);

    [[nodiscard]] const std::vector<Shell>& shells() const noexcept { return shells_; }
    [[nodiscard]] std::size_t function_count() const noexcept { return function_count_; }

private:
    std::vector<Shell> shells_;
    std::size_t function_count_ = 0;
};

} // namespace chem::kumo
