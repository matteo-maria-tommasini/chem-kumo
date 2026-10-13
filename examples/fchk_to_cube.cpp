#include <chem/kumo.hpp>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=3) { std::cerr << "Usage: kumo-example input.fchk output.cube\n"; return 2; }
    try {
        auto molecule=chem::kumo::Molecule::from_fchk(argv[1]);
        auto grid=chem::kumo::Grid::around(molecule.atoms(),4.0,0.25);
        molecule.density(grid).write_cube(argv[2]);
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
