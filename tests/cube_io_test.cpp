#include "chem/kumo/cube.hpp"
#include "chem/kumo/field_sampler.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition) {
    if (!condition) throw std::runtime_error("Cube I/O regression failed");
}
void expect_failure(const std::filesystem::path& path) {
    bool failed = false;
    try { (void)chem::kumo::read_cube(path); }
    catch (const std::exception&) { failed = true; }
    require(failed);
}
}

int main() {
    namespace fs = std::filesystem;
    const fs::path base = fs::temp_directory_path() / "kumo-cube-io-regression.cube";
    const fs::path bad = fs::temp_directory_path() / "kumo-cube-io-invalid.cube";
    try {
        chem::kumo::Grid grid;
        grid.counts = Eigen::Vector3i(2, 3, 4);
        grid.origin = Eigen::Vector3d(0.2, -0.3, 0.5);
        grid.steps << 0.2, 0.01, 0.03,
                      0.0, 0.3, 0.02,
                      0.0, 0.0, 0.4;
        std::vector<double> samples(24);
        for (std::size_t i = 0; i < samples.size(); ++i)
            samples[i] = static_cast<double>(i) / 7.0;
        chem::kumo::CubeData cube{
            chem::kumo::Field(grid, {{26, 16.0, Eigen::Vector3d(0.1, 0.2, 0.3)}}, samples),
            {"first comment", "second comment", true}};
        chem::kumo::write_cube(cube, base);
        auto loaded = chem::kumo::read_cube(base);
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
        const auto point = grid.origin + grid.steps * Eigen::Vector3d(0.4, 1.2, 2.3);
        chem::kumo::FieldSampler sampler(loaded.field);
        require(sampler.sample(point).inside);

        { std::ofstream out(bad); out << "a\nb\n0 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n1 2 3\n"; }
        expect_failure(bad);
        { std::ofstream out(bad); out << "a\nb\n0 0 0 0 2\n2 1 0 0\n2 0 1 0\n2 0 0 1\n"; }
        expect_failure(bad);
        { std::ofstream out(bad); out << "a\nb\n-1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n"; }
        expect_failure(bad);
    } catch (...) {
        fs::remove(base); fs::remove(bad); throw;
    }
    fs::remove(base); fs::remove(bad);
}
