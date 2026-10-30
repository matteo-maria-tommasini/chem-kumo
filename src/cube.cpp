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
#include <array>
#include <memory>
#include <streambuf>
#include <zlib.h>

namespace chem::kumo {
namespace {
constexpr double angstrom_to_bohr = 1.8897261246257702;

std::istringstream header_row(std::istream& input) {
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
// A bounded streaming gzip decoder. Errors are retained separately because
// formatted iostream extraction may translate streambuf exceptions into EOF.
class GzipBuffer final : public std::streambuf {
public:
    explicit GzipBuffer(const std::filesystem::path& path) : source_(path, std::ios::binary) {
        if (!source_) throw std::runtime_error("cannot open Cube input: " + path.string());
        if (inflateInit2(&z_, MAX_WBITS + 16) != Z_OK)
            throw std::runtime_error("cannot initialize Cube gzip decoder");
        initialized_ = true;
        setg(output_.data(), output_.data(), output_.data());
    }
    ~GzipBuffer() override { if (initialized_) inflateEnd(&z_); }
    GzipBuffer(const GzipBuffer&) = delete;
    GzipBuffer& operator=(const GzipBuffer&) = delete;
    void validate() const {
        if (!error_.empty()) throw std::runtime_error(error_);
        if (!finished_) throw std::runtime_error("truncated gzip Cube stream");
    }
protected:
    int_type underflow() override {
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
        if (finished_ || !error_.empty()) return traits_type::eof();
        for (;;) {
            if (!z_.avail_in) {
                source_.read(reinterpret_cast<char*>(input_.data()), input_.size());
                const auto n = source_.gcount();
                if (n == 0) { error_ = "truncated gzip Cube stream"; return traits_type::eof(); }
                z_.next_in = input_.data();
                z_.avail_in = static_cast<uInt>(n);
            }
            z_.next_out = reinterpret_cast<Bytef*>(output_.data());
            z_.avail_out = static_cast<uInt>(output_.size());
            const int status = inflate(&z_, Z_NO_FLUSH);
            const auto produced = output_.size() - z_.avail_out;
            if (status == Z_STREAM_END) {
                finished_ = true;
                if (z_.avail_in || source_.peek() != std::char_traits<char>::eof())
                    error_ = "trailing data after gzip Cube member";
            } else if (status != Z_OK) {
                error_ = "invalid gzip Cube stream or checksum";
            }
            if (produced) {
                setg(output_.data(), output_.data(), output_.data() + produced);
                return traits_type::to_int_type(*gptr());
            }
            if (finished_ || !error_.empty()) return traits_type::eof();
        }
    }
private:
    std::ifstream source_;
    z_stream z_{};
    bool initialized_ = false;
    bool finished_ = false;
    std::string error_;
    std::array<unsigned char, 65536> input_{};
    std::array<char, 65536> output_{};
};

CubeData parse_cube(std::istream& input);
} // namespace

CubeData read_cube(const std::filesystem::path& path) {
    std::ifstream probe(path, std::ios::binary);
    if (!probe) throw std::runtime_error("cannot open Cube input: " + path.string());
    unsigned char signature[2]{};
    probe.read(reinterpret_cast<char*>(signature), 2);
    const bool gzip = probe.gcount() == 2 && signature[0] == 0x1f && signature[1] == 0x8b;
    probe.clear();
    probe.seekg(0);
    if (!gzip) return parse_cube(probe);
    probe.close();
    GzipBuffer buffer(path);
    std::istream decoded(&buffer);
    CubeData result = parse_cube(decoded);
    // Consume to the physical gzip end to verify trailer CRC32/ISIZE.
    decoded.clear();
    std::array<char, 4096> discard{};
    while (decoded.read(discard.data(), discard.size()) || decoded.gcount()) {}
    buffer.validate();
    return result;
}

namespace {
CubeData parse_cube(std::istream& input) {
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
} // namespace

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
