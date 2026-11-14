#include <chem/kumo.hpp>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: kumo-example input.fchk output.cube\n";
        return 2;
    }
    try {
        auto molecule = chem::kumo::Molecule::from_fchk(argv[1]);
        auto grid = chem::kumo::Grid::around(molecule.atoms(), 4.0, 0.25);
        const auto density = chem::kumo::density::total.value();
        const auto results =
            chem::kumo::Evaluator(molecule).evaluate(grid, chem::kumo::Request().add(density));
        const auto values = results.get(density).component();
        chem::kumo::Field(grid, molecule.atoms(), {values.begin(), values.end()})
            .write_cube(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
