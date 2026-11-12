#include "chem/kumo/cube.hpp"
#include "chem/kumo/field_sampler.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>
#include <iterator>

namespace {
void require(bool condition) {
    if (!condition) throw std::runtime_error("Cube I/O regression failed");
}

// Reading `path` must fail; for gzip damage, the message must name the gzip
// problem rather than a Cube format error caused by the missing data.
void expect_failure(const std::filesystem::path& path, const std::string& expected = "") {
    try {
        (void)chem::kumo::read_cube(path);
    } catch (const std::exception& e) {
        if (std::string(e.what()).find(expected) == std::string::npos)
            throw std::runtime_error("expected an error containing '" + expected + "' for " +
                                     path.string() + ", got: " + e.what());
        return;
    }
    throw std::runtime_error("invalid Cube accepted: " + path.string());
}
} // namespace

int main() {
    namespace fs = std::filesystem;
    const fs::path base = fs::temp_directory_path() / "kumo-cube-io-regression.cube";
    const fs::path zipped = base.string() + ".gz";
    const fs::path corrupt = base.string() + ".corrupt.gz";
    const fs::path bad = fs::temp_directory_path() / "kumo-cube-io-invalid.cube";
    try {
        chem::kumo::Grid grid;
        grid.counts = Eigen::Vector3i(2, 3, 4);
        grid.origin = Eigen::Vector3d(0.2, -0.3, 0.5);
        grid.steps << 0.2, 0.01, 0.03, 0.0, 0.3, 0.02, 0.0, 0.0, 0.4;
        std::vector<double> samples(24);
        for (std::size_t i = 0; i < samples.size(); ++i)
            samples[i] = static_cast<double>(i) / 7.0;
        chem::kumo::CubeData cube{
            chem::kumo::Field(grid, {{26, 16.0, Eigen::Vector3d(0.1, 0.2, 0.3)}}, samples),
            {"first comment", "second comment", true}};
        chem::kumo::write_cube(cube, base);
        // Both writers reject nonfinite samples (fail closed): Cube readers,
        // read_cube() included, would not accept them.
        for (const double invalid : {std::nan(""), HUGE_VAL, -HUGE_VAL}) {
            std::vector<double> broken = samples;
            broken[5] = invalid;
            const chem::kumo::Field field(grid, cube.field.atoms(), broken);
            bool rejected = false;
            try {
                field.write_cube(bad);
            } catch (const std::invalid_argument&) {
                rejected = true;
            }
            if (!rejected)
                throw std::runtime_error("Field::write_cube accepted a nonfinite sample");
        }
        // Compress the exact serialized Cube; gzip and plain must agree.
        {
            std::ifstream in(base, std::ios::binary);
            gzFile out = gzopen(zipped.string().c_str(), "wb");
            require(out != nullptr);
            char chunk[4096];
            while (in.read(chunk, sizeof(chunk)) || in.gcount())
                require(gzwrite(out, chunk, static_cast<unsigned>(in.gcount())) == in.gcount());
            require(gzclose(out) == Z_OK);
        }
        auto loaded = chem::kumo::read_cube(base);
        auto gzloaded = chem::kumo::read_cube(zipped);
        require(gzloaded.metadata.first_comment == loaded.metadata.first_comment);
        require(gzloaded.metadata.second_comment == loaded.metadata.second_comment);
        require(gzloaded.metadata.coordinates_in_angstrom ==
                loaded.metadata.coordinates_in_angstrom);
        require((gzloaded.field.grid().origin - loaded.field.grid().origin).norm() == 0);
        require((gzloaded.field.grid().steps - loaded.field.grid().steps).norm() == 0);
        require(gzloaded.field.grid().counts == loaded.field.grid().counts);
        require(gzloaded.field.atoms().size() == loaded.field.atoms().size());
        require(gzloaded.field.values().size() == loaded.field.values().size());
        for (std::size_t i = 0; i < loaded.field.values().size(); ++i)
            require(gzloaded.field.values()[i] == loaded.field.values()[i]);
        {
            std::ifstream in(zipped, std::ios::binary);
            std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
            require(bytes.size() > 16);
            bytes[bytes.size() - 5] ^= 0x40; // CRC32 trailer corruption.
            std::ofstream out(corrupt, std::ios::binary);
            out.write(bytes.data(), bytes.size());
        }
        expect_failure(corrupt, "invalid gzip data or checksum");
        {
            std::ifstream in(zipped, std::ios::binary);
            std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
            std::ofstream out(corrupt, std::ios::binary | std::ios::trunc);
            out.write(bytes.data(), bytes.size() - 6);
        }
        expect_failure(corrupt, "truncated gzip");
        {
            // Truncation inside the compressed data: the parser sees a short
            // Cube, but the reported cause must be the truncated gzip stream.
            std::ifstream in(zipped, std::ios::binary);
            std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
            std::ofstream out(corrupt, std::ios::binary | std::ios::trunc);
            out.write(bytes.data(), static_cast<std::streamsize>(bytes.size() / 2));
        }
        expect_failure(corrupt, "truncated gzip");
        {
            std::ifstream in(zipped, std::ios::binary);
            std::ofstream out(corrupt, std::ios::binary | std::ios::trunc);
            out << in.rdbuf() << "garbage";
        }
        expect_failure(corrupt, "trailing data after gzip member");
        require(loaded.metadata.first_comment == "first comment");
        require(loaded.metadata.second_comment == "second comment");
        require(loaded.metadata.coordinates_in_angstrom);
        require(loaded.field.atoms().size() == 1);
        require(loaded.field.atoms()[0].atomic_number == 26);
        require(loaded.field.atoms()[0].nuclear_charge == 16.0);
        require((loaded.field.grid().steps - grid.steps).norm() < 1e-11);
        require((loaded.field.grid().origin - grid.origin).norm() < 1e-11);
        require((loaded.field.atoms()[0].position - cube.field.atoms()[0].position).norm() < 1e-11);
        require(loaded.field.values().size() == samples.size());
        for (std::size_t i = 0; i < samples.size(); ++i)
            require(std::abs(loaded.field.values()[i] - samples[i]) < 1e-11);
        // A concrete Vector3d, not auto: an auto Eigen expression would keep a
        // reference to the temporary Vector3d after this statement.
        const Eigen::Vector3d point = grid.origin + grid.steps * Eigen::Vector3d(0.4, 1.2, 2.3);
        chem::kumo::FieldSampler sampler(loaded.field);
        require(sampler.sample(point).inside);

        {
            std::ofstream out(bad);
            out << "a\nb\n0 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n1 2 3\n";
        }
        expect_failure(bad);
        {
            std::ofstream out(bad);
            out << "a\nb\n0 0 0 0 2\n2 1 0 0\n2 0 1 0\n2 0 0 1\n";
        }
        expect_failure(bad);
        {
            std::ofstream out(bad);
            out << "a\nb\n-1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n";
        }
        expect_failure(bad);
    } catch (...) {
        fs::remove(base);
        fs::remove(bad);
        fs::remove(zipped);
        fs::remove(corrupt);
        throw;
    }
    fs::remove(base);
    fs::remove(bad);
    fs::remove(zipped);
    fs::remove(corrupt);
}
