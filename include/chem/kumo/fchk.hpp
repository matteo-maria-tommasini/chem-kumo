#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace chem::kumo {

/// Typed reader for named records in a Gaussian formatted checkpoint file.
/// Missing records or mismatched types are errors; no binary CHK support.
class FormattedCheckpoint {
public:
    explicit FormattedCheckpoint(const std::filesystem::path& path);

    [[nodiscard]] bool contains(const std::string& label) const;
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
