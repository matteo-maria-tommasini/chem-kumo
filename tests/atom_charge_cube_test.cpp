#include "chem/kumo/density.hpp"
#include "test_support.hpp"

#include <Eigen/Core>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void verify_charge(const std::filesystem::path& path, int expected_z, double expected_charge,
                   bool orbital) {
    const std::string name = path.filename().string() + ": ";
    std::ifstream input(path);
    if (!kumo_test::check(static_cast<bool>(input), name + "the written Cube file can be opened"))
        return;
    std::string line;
    std::getline(input, line);
    std::getline(input, line);
    int atom_count = 0;
    double x = 0.0, y = 0.0, z = 0.0;
    if (!kumo_test::check(static_cast<bool>(input >> atom_count >> x >> y >> z),
                          name + "readable atom count and origin line"))
        return;
    const int expected_count = orbital ? -1 : 1;
    kumo_test::check(atom_count == expected_count, name + "atom count field is " +
                                                       std::to_string(expected_count) + " (got " +
                                                       std::to_string(atom_count) + ")");
    for (int axis = 0; axis < 3; ++axis) {
        int count = 0;
        if (!kumo_test::check(static_cast<bool>(input >> count >> x >> y >> z),
                              name + "readable grid line " + std::to_string(axis)))
            return;
        kumo_test::check(count == 2, name + "grid axis " + std::to_string(axis) +
                                         " has 2 points (got " + std::to_string(count) + ")");
    }
    int atomic_number = 0;
    double nuclear_charge = 0.0;
    if (!kumo_test::check(
            static_cast<bool>(input >> atomic_number >> nuclear_charge >> x >> y >> z),
            name + "readable atom line"))
        return;
    kumo_test::check(atomic_number == expected_z, name + "atomic number " +
                                                      std::to_string(expected_z) + " (got " +
                                                      std::to_string(atomic_number) + ")");
    kumo_test::check(std::abs(nuclear_charge - expected_charge) < 1e-6,
                     name + "nuclear charge " + kumo_test::num(expected_charge) +
                         " within 1e-6 (got " + kumo_test::num(nuclear_charge, 17) + ")");
}

} // namespace

int main() {
    const auto density_path =
        std::filesystem::temp_directory_path() / "kumo-charge-density-test.cube";
    const auto orbital_path =
        std::filesystem::temp_directory_path() / "kumo-charge-orbital-test.cube";
    const int result = kumo_test::run([&] {
        chem::kumo::Grid grid;
        grid.counts = Eigen::Vector3i(2, 2, 2);
        const std::vector<chem::kumo::Atom> atoms{{26, 16.0, Eigen::Vector3d(0.0, 0.0, 0.0)}};
        chem::kumo::Field field(grid, atoms, std::vector<double>(8, 0.5));
        kumo_test::check(field.atoms()[0].atomic_number == 26, "Field keeps the atomic number 26");
        kumo_test::check(field.atoms()[0].nuclear_charge == 16.0,
                         "Field keeps the nuclear charge 16 (an ECP core charge)");
        field.write_cube(density_path);
        field.write_orbital_cube(orbital_path, 1);
        verify_charge(density_path, 26, 16.0, false);
        verify_charge(orbital_path, 26, 16.0, true);
    });
    std::filesystem::remove(density_path);
    std::filesystem::remove(orbital_path);
    return result;
}
