#pragma once

#include "chem/kumo/density.hpp"

#include <filesystem>
#include <string>

namespace chem::kumo {

/// Header information of a Cube file that is not part of the numerical
/// field: the two comment lines and the length unit of the header.
struct CubeMetadata {
    std::string first_comment;  ///< first header line
    std::string second_comment; ///< second header line
    /// True when the header used angstrom (negative grid counts); the Field
    /// is in bohr either way.
    bool coordinates_in_angstrom = false;
};

/// One scalar Cube dataset: the field (always in bohr) and its header
/// information.
struct CubeData {
    Field field;
    CubeMetadata metadata;
};

/// Read a scalar Gaussian Cube file, plain or gzip-compressed (detected from
/// its content). The atom count must not be negative (orbital Cubes are
/// rejected), and the optional number of values per point (NVAL) after the
/// origin must be 1. Angstrom headers are converted to bohr. Throws
/// std::runtime_error if the file cannot be read or is malformed.
[[nodiscard]] CubeData read_cube(const std::filesystem::path& path);

/// Write a scalar Cube file with the comment lines and the length unit of
/// `cube.metadata`, values with 13 significant digits. Throws
/// std::invalid_argument for a comment with a line break, a nonfinite value
/// or an invalid atom, and std::runtime_error if the file cannot be written.
void write_cube(const CubeData& cube, const std::filesystem::path& path);

} // namespace chem::kumo
