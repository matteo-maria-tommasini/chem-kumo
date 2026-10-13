#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace chem::kumo::detail {

inline std::size_t checked_add(std::size_t a, std::size_t b, const char* what) {
    if (a > std::numeric_limits<std::size_t>::max() - b) {
        throw std::overflow_error(std::string(what) + " overflows size_t");
    }
    return a + b;
}

inline std::size_t checked_mul(std::size_t a, std::size_t b, const char* what) {
    if (a != 0U && b > std::numeric_limits<std::size_t>::max() / a) {
        throw std::overflow_error(std::string(what) + " overflows size_t");
    }
    return a * b;
}

inline std::size_t checked_size(std::int64_t value, const char* what) {
    if (value < 0) {
        throw std::runtime_error(std::string(what) + " is negative");
    }
    const auto unsigned_value = static_cast<std::uint64_t>(value);
    if (unsigned_value > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error(std::string(what) + " does not fit size_t");
    }
    return static_cast<std::size_t>(unsigned_value);
}

inline std::size_t checked_positive_size(std::int64_t value, const char* what) {
    const auto result = checked_size(value, what);
    if (result == 0U) {
        throw std::runtime_error(std::string(what) + " must be positive");
    }
    return result;
}

} // namespace chem::kumo::detail
