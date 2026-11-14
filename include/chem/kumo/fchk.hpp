#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace chem::kumo {

/// Reader of the named records (scalars and arrays, integer or real) of a
/// Gaussian formatted checkpoint file (FCHK). Asking for a missing record or
/// for a record of another type is an error. Binary checkpoint files (.chk)
/// are not supported.
class FormattedCheckpoint {
public:
    /// Read and parse a plain or gzip-compressed FCHK (detected from its
    /// content). Throws std::runtime_error if the file cannot be read or a
    /// record is malformed or truncated.
    explicit FormattedCheckpoint(const std::filesystem::path& path);

    /// Whether a record with this exact label exists.
    [[nodiscard]] bool contains(const std::string& label) const;
    /// Value of an integer scalar, real scalar, integer array or real array
    /// record. Each throws std::runtime_error if the record is missing or
    /// has another type (scalar or array, integer or real).
    [[nodiscard]] long long integer(const std::string& label) const;
    [[nodiscard]] double real(const std::string& label) const;
    [[nodiscard]] const std::vector<long long>& integers(const std::string& label) const;
    [[nodiscard]] const std::vector<double>& reals(const std::string& label) const;

private:
    enum class Kind { IntegerScalar, RealScalar, IntegerArray, RealArray };

    struct Record {
        Kind kind{};
        long long integer_scalar = 0;
        double real_scalar = 0.0;
        std::vector<long long> integer_array;
        std::vector<double> real_array;
    };

    [[nodiscard]] const Record& record(const std::string& label) const;

    std::unordered_map<std::string, Record> records_;
};

} // namespace chem::kumo
