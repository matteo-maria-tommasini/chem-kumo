#include "chem/kumo/cube.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chem::kumo {
namespace {
constexpr double angstrom_to_bohr = 1.8897261246257702;

std::istringstream header_row(std::ifstream& input) {
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("truncated Cube header");
    std::replace(line.begin(), line.end(), 'D', 'E');
    std::replace(line.begin(), line.end(), 'd', 'e');
    return std::istringstream(line);
}

void end_row(std::istringstream& row) {
    row >> std::ws;
    if (!row.eof()) throw std::runtime_error("malformed Cube header");
}

std::size_t sample_count(const Grid& grid) {
    std::size_t count = 1;
    for (int axis = 0; axis < 3; ++axis) {
        const int n = grid.counts[axis];
        if (n < 2 || static_cast<std::size_t>(n) > 1000000000ULL / count)
            throw std::runtime_error("invalid Cube dimensions");
        count *= static_cast<std::size_t>(n);
    }
    return count;
}

void check_comments(const CubeMetadata& metadata) {
    if (metadata.first_comment.find_first_of("\r\n") != std::string::npos ||
        metadata.second_comment.find_first_of("\r\n") != std::string::npos)
        throw std::invalid_argument("Cube comments must be single lines");
}
} // namespace

CubeData read_cube(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Cube input: " + path.string());
    CubeMetadata metadata;
    if (!std::getline(input, metadata.first_comment) ||
        !std::getline(input, metadata.second_comment))
        throw std::runtime_error("missing Cube comments");

    Grid grid;
    auto header = header_row(input);
    long long atom_count = 0;
    if (!(header >> atom_count >> grid.origin[0] >> grid.origin[1] >> grid.origin[2]) ||
        atom_count < 0 || atom_count > 1000000)
        throw std::runtime_error("invalid Cube atom count or unsupported orbital Cube");
    header >> std::ws;
    if (!header.eof()) {
        int nval = 0;
        if (!(header >> nval) || nval != 1)
            throw std::runtime_error("only scalar Cube NVAL=1 is supported");
    }
    end_row(header);

    for (int axis = 0; axis < 3; ++axis) {
        auto row = header_row(input);
        long long signed_count = 0;
        Eigen::Vector3d step;
        if (!(row >> signed_count >> step[0] >> step[1] >> step[2]) ||
            signed_count < -1000000000LL || signed_count > 1000000000LL ||
            (signed_count >= -1 && signed_count <= 1))
            throw std::runtime_error("invalid Cube grid axis");
        end_row(row);
        if (axis == 0) metadata.coordinates_in_angstrom = signed_count < 0;
        if ((signed_count < 0) != metadata.coordinates_in_angstrom)
            throw std::runtime_error("mixed Cube coordinate units");
        grid.counts[axis] = static_cast<int>(std::abs(signed_count));
        grid.steps.col(axis) = step;
    }
    const std::size_t count = sample_count(grid);

    std::vector<Atom> atoms;
    atoms.reserve(static_cast<std::size_t>(atom_count));
    for (long long i = 0; i < atom_count; ++i) {
        auto row = header_row(input);
        Atom atom;
        if (!(row >> atom.atomic_number >> atom.nuclear_charge >>
              atom.position[0] >> atom.position[1] >> atom.position[2]))
            throw std::runtime_error("invalid Cube atom record");
        end_row(row);
        if (atom.atomic_number < 0 || atom.atomic_number > 118 ||
            !std::isfinite(atom.nuclear_charge) || !atom.position.allFinite())
            throw std::runtime_error("invalid Cube atom");
        atoms.push_back(atom);
    }
    if (metadata.coordinates_in_angstrom) {
        grid.origin *= angstrom_to_bohr;
        grid.steps *= angstrom_to_bohr;
        for (auto& atom : atoms) atom.position *= angstrom_to_bohr;
    }

    std::vector<double> values;
    values.reserve(count);
    std::string token;
    while (input >> token) {
        if (values.size() == count) throw std::runtime_error("extra Cube samples");
        std::replace(token.begin(), token.end(), 'D', 'E');
        std::replace(token.begin(), token.end(), 'd', 'e');
        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (end != token.c_str() + token.size() || !std::isfinite(value))
            throw std::runtime_error("invalid Cube sample");
        values.push_back(value);
    }
    if (input.bad()) throw std::runtime_error("Cube input read failure");
    if (values.size() != count) throw std::runtime_error("truncated Cube samples");
    return {Field(std::move(grid), std::move(atoms), std::move(values)),
            std::move(metadata)};
}

void write_cube(const CubeData& cube, const std::filesystem::path& path) {
    check_comments(cube.metadata);
    const auto& grid = cube.field.grid();
    const auto& atoms = cube.field.atoms();
    const auto values = cube.field.values();
    if (values.size() != sample_count(grid))
        throw std::invalid_argument("inconsistent Cube field size");
    const double factor = cube.metadata.coordinates_in_angstrom ? 1.0 / angstrom_to_bohr : 1.0;
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open Cube output: " + path.string());
    output << cube.metadata.first_comment << '\n' << cube.metadata.second_comment << '\n';
    output << std::scientific << std::setprecision(12);
    output << atoms.size();
    for (int k = 0; k < 3; ++k) output << ' ' << grid.origin[k] * factor;
    output << '\n';
    for (int axis = 0; axis < 3; ++axis) {
        output << (cube.metadata.coordinates_in_angstrom ? -grid.counts[axis] : grid.counts[axis]);
        for (int k = 0; k < 3; ++k) output << ' ' << grid.steps(k, axis) * factor;
        output << '\n';
    }
    for (const auto& atom : atoms) {
        if (atom.atomic_number < 0 || atom.atomic_number > 118 ||
            !std::isfinite(atom.nuclear_charge) || !atom.position.allFinite())
            throw std::invalid_argument("invalid Cube atom");
        output << atom.atomic_number << ' ' << atom.nuclear_charge;
        for (int k = 0; k < 3; ++k) output << ' ' << atom.position[k] * factor;
        output << '\n';
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) throw std::invalid_argument("nonfinite Cube sample");
        output << values[i] << ((i + 1) % 6 == 0 ? '\n' : ' ');
    }
    if (values.size() % 6 != 0) output << '\n';
    output.flush();
    if (!output) throw std::runtime_error("Cube output write failed");
}
} // namespace chem::kumo
