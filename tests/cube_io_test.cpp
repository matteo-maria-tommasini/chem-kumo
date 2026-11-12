#include "chem/kumo/cube.hpp"
#include "chem/kumo/field_sampler.hpp"
#include "test_support.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>
#include <iterator>

namespace {
// Reading `path` must fail; for gzip damage, the message must name the gzip
// problem rather than a Cube format error caused by the missing data.
void expect_failure(const std::filesystem::path& path, const std::string& expected = "") {
    try {
        (void)chem::kumo::read_cube(path);
    } catch (const std::exception& e) {
        kumo_test::check(std::string(e.what()).find(expected) != std::string::npos,
                         "expected an error containing '" + expected + "' for " + path.string() +
                             ", got: " + e.what());
        return;
    }
    kumo_test::check(false, "read_cube rejects the invalid Cube " + path.string() +
                                " (no exception thrown)");
}
} // namespace

int main() {
    namespace fs = std::filesystem;
    const fs::path base = fs::temp_directory_path() / "kumo-cube-io-regression.cube";
    const fs::path zipped = base.string() + ".gz";
    const fs::path corrupt = base.string() + ".corrupt.gz";
    const fs::path bad = fs::temp_directory_path() / "kumo-cube-io-invalid.cube";
    const int result = kumo_test::run([&] {
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
            kumo_test::check(rejected, "Field::write_cube rejects the nonfinite sample " +
                                           kumo_test::num(invalid));
        }
        // Compress the exact serialized Cube; gzip and plain must agree.
        {
            std::ifstream in(base, std::ios::binary);
            gzFile out = gzopen(zipped.string().c_str(), "wb");
            if (out == nullptr) throw std::runtime_error("cannot create " + zipped.string());
            char chunk[4096];
            while (in.read(chunk, sizeof(chunk)) || in.gcount())
                if (gzwrite(out, chunk, static_cast<unsigned>(in.gcount())) != in.gcount())
                    throw std::runtime_error("cannot write " + zipped.string());
            if (gzclose(out) != Z_OK) throw std::runtime_error("cannot close " + zipped.string());
        }
        auto loaded = chem::kumo::read_cube(base);
        auto gzloaded = chem::kumo::read_cube(zipped);
        const std::string gz = "gzip and plain Cube readers agree: ";
        kumo_test::check(gzloaded.metadata.first_comment == loaded.metadata.first_comment,
                         gz + "first comment");
        kumo_test::check(gzloaded.metadata.second_comment == loaded.metadata.second_comment,
                         gz + "second comment");
        kumo_test::check(gzloaded.metadata.coordinates_in_angstrom ==
                             loaded.metadata.coordinates_in_angstrom,
                         gz + "coordinate units");
        kumo_test::check((gzloaded.field.grid().origin - loaded.field.grid().origin).norm() == 0,
                         gz + "origin, exactly");
        kumo_test::check((gzloaded.field.grid().steps - loaded.field.grid().steps).norm() == 0,
                         gz + "steps, exactly");
        kumo_test::check(gzloaded.field.grid().counts == loaded.field.grid().counts,
                         gz + "grid counts");
        kumo_test::check(gzloaded.field.atoms().size() == loaded.field.atoms().size(),
                         gz + "atom count");
        if (kumo_test::check(gzloaded.field.values().size() == loaded.field.values().size(),
                             gz + "sample count"))
            for (std::size_t i = 0; i < loaded.field.values().size(); ++i)
                kumo_test::check(gzloaded.field.values()[i] == loaded.field.values()[i],
                                 gz + "sample " + std::to_string(i) + ", exactly");
        {
            std::ifstream in(zipped, std::ios::binary);
            std::vector<char> bytes((std::istreambuf_iterator<char>(in)), {});
            if (bytes.size() <= 16)
                throw std::runtime_error("unexpectedly short gzip file " + zipped.string());
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
        const std::string rt = "Cube write/read round trip: ";
        kumo_test::check(loaded.metadata.first_comment == "first comment",
                         rt + "first comment is 'first comment' (got '" +
                             loaded.metadata.first_comment + "')");
        kumo_test::check(loaded.metadata.second_comment == "second comment",
                         rt + "second comment is 'second comment' (got '" +
                             loaded.metadata.second_comment + "')");
        kumo_test::check(loaded.metadata.coordinates_in_angstrom,
                         rt + "coordinates stay in angstrom");
        if (kumo_test::check(loaded.field.atoms().size() == 1,
                             rt + "one atom (got " + std::to_string(loaded.field.atoms().size()) +
                                 ")")) {
            const auto& atom = loaded.field.atoms()[0];
            kumo_test::check(atom.atomic_number == 26, rt + "atomic number 26 (got " +
                                                           std::to_string(atom.atomic_number) +
                                                           ")");
            kumo_test::check(atom.nuclear_charge == 16.0, rt + "nuclear charge 16 (got " +
                                                              kumo_test::num(atom.nuclear_charge) +
                                                              ")");
            const double position_error = (atom.position - cube.field.atoms()[0].position).norm();
            kumo_test::check(position_error < 1e-11, rt + "atom position within 1e-11 (error " +
                                                         kumo_test::num(position_error) + ")");
        }
        const double steps_error = (loaded.field.grid().steps - grid.steps).norm();
        kumo_test::check(steps_error < 1e-11,
                         rt + "steps within 1e-11 (error " + kumo_test::num(steps_error) + ")");
        const double origin_error = (loaded.field.grid().origin - grid.origin).norm();
        kumo_test::check(origin_error < 1e-11,
                         rt + "origin within 1e-11 (error " + kumo_test::num(origin_error) + ")");
        if (kumo_test::check(loaded.field.values().size() == samples.size(),
                             rt + "sample count " + std::to_string(samples.size())))
            for (std::size_t i = 0; i < samples.size(); ++i)
                kumo_test::check(std::abs(loaded.field.values()[i] - samples[i]) < 1e-11,
                                 rt + "sample " + std::to_string(i) + " within 1e-11 (expected " +
                                     kumo_test::num(samples[i], 17) + ", got " +
                                     kumo_test::num(loaded.field.values()[i], 17) + ")");
        // A concrete Vector3d, not auto: an auto Eigen expression would keep a
        // reference to the temporary Vector3d after this statement.
        const Eigen::Vector3d point = grid.origin + grid.steps * Eigen::Vector3d(0.4, 1.2, 2.3);
        chem::kumo::FieldSampler sampler(loaded.field);
        kumo_test::check(sampler.sample(point).inside,
                         "a point inside the read grid is sampled as inside");

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
    });
    fs::remove(base);
    fs::remove(bad);
    fs::remove(zipped);
    fs::remove(corrupt);
    return result;
}
