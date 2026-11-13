#include "chem/kumo/cube.hpp"
#include "test_support.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

// Absolute tolerance; `expected` is always a finite literal here.
void near(double actual, double expected, double tolerance, const std::string& message) {
    kumo_test::check_close(actual, expected, 0.0, tolerance, message);
}

struct Reference {
    const char* filename;
    const char* first_comment;
    const char* second_comment;
    std::array<int, 3> counts;
    std::size_t atoms;
    int first_atomic_number;
    double first_charge;
    std::array<double, 3> first_position;
    std::array<double, 3> origin;
    std::array<double, 3> sample_values;
};

void check_reference(const std::filesystem::path& directory, const Reference& reference) {
    const auto cube = chem::kumo::read_cube(directory / reference.filename);
    const std::string file = std::string(reference.filename) + ": ";
    const auto& grid = cube.field.grid();
    const auto& atoms = cube.field.atoms();
    const auto values = cube.field.values();

    kumo_test::check(cube.metadata.first_comment == reference.first_comment,
                     file + "first comment is '" + reference.first_comment + "' (got '" +
                         cube.metadata.first_comment + "')");
    kumo_test::check(cube.metadata.second_comment == reference.second_comment,
                     file + "second comment is '" + reference.second_comment + "' (got '" +
                         cube.metadata.second_comment + "')");
    kumo_test::check(!cube.metadata.coordinates_in_angstrom, file + "Gaussian reference uses bohr");
    if (!kumo_test::check(atoms.size() == reference.atoms, file + std::to_string(reference.atoms) +
                                                               " atoms (got " +
                                                               std::to_string(atoms.size()) + ")"))
        return;
    for (int axis = 0; axis < 3; ++axis) {
        const std::string a = " (axis " + std::to_string(axis) + ")";
        kumo_test::check(grid.counts[axis] == reference.counts[axis],
                         file + "grid count " + std::to_string(reference.counts[axis]) + a +
                             " (got " + std::to_string(grid.counts[axis]) + ")");
        near(grid.origin[axis], reference.origin[axis], 1e-12, file + "grid origin" + a);
        for (int component = 0; component < 3; ++component)
            near(grid.steps(component, axis), axis == component ? 0.166667 : 0.0, 1e-12,
                 file + "grid step component " + std::to_string(component) + a);
        near(atoms.front().position[axis], reference.first_position[axis], 1e-12,
             file + "first atom position" + a);
    }
    kumo_test::check(atoms.front().atomic_number == reference.first_atomic_number,
                     file + "first atomic number " + std::to_string(reference.first_atomic_number) +
                         " (got " + std::to_string(atoms.front().atomic_number) + ")");
    near(atoms.front().nuclear_charge, reference.first_charge, 1e-12,
         file + "first nuclear charge");

    const std::size_t expected_count = static_cast<std::size_t>(reference.counts[0]) *
                                       static_cast<std::size_t>(reference.counts[1]) *
                                       static_cast<std::size_t>(reference.counts[2]);
    if (!kumo_test::check(values.size() == expected_count, file + std::to_string(expected_count) +
                                                               " samples (got " +
                                                               std::to_string(values.size()) + ")"))
        return;
    const std::array<std::size_t, 3> indices{0, expected_count / 2, expected_count - 1};
    for (std::size_t i = 0; i < indices.size(); ++i) {
        const double expected = reference.sample_values[i];
        near(values[indices[i]], expected, 1e-12 + 1e-10 * std::abs(expected),
             file + "Gaussian reference sample " + std::to_string(indices[i]));
    }
}

} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, 1, "DATA_DIRECTORY")) return 2;
    return kumo_test::run([&] {
        // Literal reference values come from independently supplied Gaussian Cube files,
        // not from chem-kumo's Cube writer or a second parser in this test.
        const std::array<Reference, 3> references{{{"o2_density.cube",
                                                    "O2 ground-state SCF density reference",
                                                    " Electron density from Total SCF Density",
                                                    {61, 61, 74},
                                                    2,
                                                    8,
                                                    8.0,
                                                    {0.0, 0.0, 1.143284},
                                                    {-4.970736, -4.970736, -6.114020},
                                                    {5.94950e-18, 0.52465, 4.44573e-18}},
                                                   {"o2_spin.cube",
                                                    "O2 ground-state SCF spin reference",
                                                    " Electron density from Spin SCF Density",
                                                    {61, 61, 74},
                                                    2,
                                                    8,
                                                    8.0,
                                                    {0.0, 0.0, 1.143284},
                                                    {-4.970736, -4.970736, -6.114020},
                                                    {1.92726e-18, -0.00166261, 1.45692e-18}},
                                                   {"acrolein_density.cube",
                                                    "Acrolein ground-state SCF density reference",
                                                    " Electron density from Total SCF Density",
                                                    {97, 102, 79},
                                                    8,
                                                    6,
                                                    6.0,
                                                    {2.272099, 2.443678, 0.0},
                                                    {-7.278292, -7.954711, -6.512752},
                                                    {3.08107e-21, 1.53016e-07, 2.85068e-17}}}};
        for (const auto& reference : references)
            check_reference(argv[1], reference);
    });
}
