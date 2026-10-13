#include "chem/kumo/basis.hpp"

#include "checked_math.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace chem::kumo {

std::size_t Shell::function_count() const {
    if (is_sp) return 4U;
    if (angular_momentum < 0) throw std::runtime_error("negative angular momentum in shell");
    const auto l = static_cast<std::size_t>(angular_momentum);
    if (representation == ShellRepresentation::Pure && angular_momentum >= 2) {
        return 2U * l + 1U;
    }
    return (l + 1U) * (l + 2U) / 2U;
}

BasisSet::BasisSet(const FormattedCheckpoint& fchk) {
    const auto& types = fchk.integers("Shell types");
    const auto& nprim = fchk.integers("Number of primitives per shell");
    const auto& centers = fchk.integers("Shell to atom map");
    const auto& exponents = fchk.reals("Primitive exponents");
    const auto& coeff = fchk.reals("Contraction coefficients");
    const auto& p_coeff = fchk.reals("P(S=P) Contraction coefficients");
    const auto& xyz = fchk.reals("Coordinates of each shell");

    if (types.size() != nprim.size() || types.size() != centers.size()) {
        throw std::runtime_error("inconsistent shell arrays in fchk");
    }
    if (xyz.size() != 3U * types.size()) {
        throw std::runtime_error("invalid shell-coordinate array in fchk");
    }
    if (exponents.size() != coeff.size() || exponents.size() != p_coeff.size()) {
        throw std::runtime_error("inconsistent primitive arrays in fchk");
    }

    const auto atom_count = detail::checked_positive_size(fchk.integer("Number of atoms"), "Number of atoms");
    const auto declared_basis = detail::checked_positive_size(fchk.integer("Number of basis functions"),
                                                      "Number of basis functions");

    shells_.reserve(types.size());
    std::size_t primitive_offset = 0;
    for (std::size_t i = 0; i < types.size(); ++i) {
        const long long raw_type = types[i];
        const bool sp = raw_type == -1;
        if (raw_type == std::numeric_limits<long long>::min()) {
            throw std::runtime_error("invalid shell type in fchk");
        }
        const long long l_raw = sp ? 1LL : (raw_type < 0 ? -raw_type : raw_type);
        if (l_raw < 0 || l_raw > static_cast<long long>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("unsupported angular momentum in fchk shell");
        }
        const int l = static_cast<int>(l_raw);
        const auto count = detail::checked_positive_size(nprim[i], "Number of primitives per shell");
        if (count > exponents.size() - std::min(primitive_offset, exponents.size())) {
            throw std::runtime_error("invalid primitive count in fchk shell");
        }
        if (centers[i] <= 0 || static_cast<std::size_t>(centers[i]) > atom_count) {
            throw std::runtime_error("invalid 1-based shell center in fchk");
        }

        Shell shell;
        shell.center = static_cast<std::size_t>(centers[i] - 1);
        shell.angular_momentum = l;
        shell.is_sp = sp;
        shell.representation = (!sp && raw_type <= -2) ? ShellRepresentation::Pure
                                                        : ShellRepresentation::Cartesian;
        shell.origin = Eigen::Vector3d(xyz[3U * i], xyz[3U * i + 1U], xyz[3U * i + 2U]);
        if (!shell.origin.allFinite()) {
            throw std::runtime_error("non-finite shell coordinate in fchk");
        }
        shell.exponents.assign(exponents.begin() + static_cast<std::ptrdiff_t>(primitive_offset),
                               exponents.begin() + static_cast<std::ptrdiff_t>(primitive_offset + count));
        shell.coefficients.assign(coeff.begin() + static_cast<std::ptrdiff_t>(primitive_offset),
                                  coeff.begin() + static_cast<std::ptrdiff_t>(primitive_offset + count));
        shell.p_coefficients.assign(p_coeff.begin() + static_cast<std::ptrdiff_t>(primitive_offset),
                                    p_coeff.begin() + static_cast<std::ptrdiff_t>(primitive_offset + count));
        for (std::size_t p = 0; p < count; ++p) {
            if (!(shell.exponents[p] > 0.0) || !std::isfinite(shell.exponents[p]) ||
                !std::isfinite(shell.coefficients[p]) || !std::isfinite(shell.p_coefficients[p])) {
                throw std::runtime_error("invalid primitive exponent/coefficient in fchk shell");
            }
        }
        primitive_offset += count;
        const auto functions = shell.function_count();
        if (functions > declared_basis - std::min(function_count_, declared_basis)) {
            throw std::runtime_error("basis function count exceeds fchk declaration");
        }
        function_count_ += functions;
        shells_.push_back(std::move(shell));
    }
    if (primitive_offset != exponents.size()) {
        throw std::runtime_error("unused primitive data after parsing fchk basis");
    }
    if (function_count_ != declared_basis) {
        throw std::runtime_error("basis size reconstructed from shells does not match fchk");
    }
}

} // namespace chem::kumo
