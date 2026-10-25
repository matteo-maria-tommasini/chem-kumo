#include "chem/kumo/density.hpp"

#include <Eigen/Core>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition) {
    if (!condition) throw std::runtime_error("atom charge Cube regression failed");
}

void verify_charge(const std::filesystem::path& path, int expected_z,
                   double expected_charge, bool orbital) {
    std::ifstream input(path);
    require(static_cast<bool>(input));
    std::string line;
    std::getline(input, line);
    std::getline(input, line);
    int atom_count = 0;
    double x = 0.0, y = 0.0, z = 0.0;
    require(static_cast<bool>(input >> atom_count >> x >> y >> z));
    require(atom_count == (orbital ? -1 : 1));
    for (int axis = 0; axis < 3; ++axis) {
        int count = 0;
        require(static_cast<bool>(input >> count >> x >> y >> z));
        require(count == 2);
    }
    int atomic_number = 0;
    double nuclear_charge = 0.0;
    require(static_cast<bool>(input >> atomic_number >> nuclear_charge >> x >> y >> z));
    require(atomic_number == expected_z);
    require(std::abs(nuclear_charge - expected_charge) < 1e-6);
}

} // namespace

int main() {
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(2, 2, 2);
    const std::vector<chem::kumo::Atom> atoms{
        {26, 16.0, Eigen::Vector3d(0.0, 0.0, 0.0)}
    };
    chem::kumo::Field field(grid, atoms, std::vector<double>(8, 0.5));
    require(field.atoms()[0].atomic_number == 26);
    require(field.atoms()[0].nuclear_charge == 16.0);
    const auto density_path = std::filesystem::temp_directory_path() / "kumo-charge-density-test.cube";
    const auto orbital_path = std::filesystem::temp_directory_path() / "kumo-charge-orbital-test.cube";
    try {
        field.write_cube(density_path);
        field.write_orbital_cube(orbital_path, 1);
        verify_charge(density_path, 26, 16.0, false);
        verify_charge(orbital_path, 26, 16.0, true);
    } catch (...) {
        std::filesystem::remove(density_path);
        std::filesystem::remove(orbital_path);
        throw;
    }
    std::filesystem::remove(density_path);
    std::filesystem::remove(orbital_path);
}
