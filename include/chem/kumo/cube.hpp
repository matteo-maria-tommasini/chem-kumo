#pragma once

#include "chem/kumo/density.hpp"

#include <filesystem>
#include <string>

namespace chem::kumo {

/// Format-specific information, separate from the numerical field.
struct CubeMetadata {
    std::string first_comment;
    std::string second_comment;
    bool coordinates_in_angstrom = false;
};

/// One scalar Cube dataset. Coordinates in Field are always in bohr.
struct CubeData {
    Field field;
    CubeMetadata metadata;
};

/// Read a scalar Cube with positive atom count and NVAL absent or equal to one.
[[nodiscard]] CubeData read_cube(const std::filesystem::path& path);

/// Write a scalar Cube, retaining comments and the original coordinate convention.
void write_cube(const CubeData& cube, const std::filesystem::path& path);

} // namespace chem::kumo
