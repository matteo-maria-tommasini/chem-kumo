// kumo: the command-line program of chem-kumo. Options are parsed in
// cli/options.cpp; each mode writes its output through one of the cli/*
// modules below, and all numerical work is done by libkumo.
#include "chem/kumo.hpp"
#include "cli/console.hpp"
#include "cli/density_output.hpp"
#include "cli/options.hpp"
#include "cli/orbital_output.hpp"
#include "cli/properties.hpp"
#include "cli/sampling_output.hpp"
#include "cli/summary.hpp"
#include "cli/surface_output.hpp"

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    using namespace kumo_cli;
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") {
            usage();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--version") {
            std::cout << "kumo " CHEM_KUMO_VERSION "\n";
            return 0;
        }
        const auto o = parse(argc, argv);
        Console console(o.verbose);
        console.say("kumo " CHEM_KUMO_VERSION);
        console.begin("Reading formatted checkpoint " + o.input.string());
        const auto mol = chem::kumo::Molecule::from_fchk(o.input);
        console.end("Read in");
        console.say("  " + molecule_summary(mol));
        console.say("  " + orbital_check_summary(mol));
        switch (o.mode) {
        case Mode::ListOrbitals: list_orbitals(mol, o); return 0;
        case Mode::SurfaceMap:
            report_surface_map(mol, o, select_density(o, mol, console), console);
            break;
        case Mode::EspSurface: report_esp_surface(mol, o, console); break;
        case Mode::AlieSurface: report_alie_surface(mol, o, console); break;
        case Mode::Points: print_points(mol, o, select_density(o, mol, console), console); break;
        case Mode::Line:
        case Mode::Plane: write_profile(mol, o, select_density(o, mol, console), console); break;
        case Mode::OrbitalCubes:
        case Mode::DensityCubes: {
            const auto grid = chem::kumo::Grid::around(mol.atoms(), o.margin, o.spacing);
            console.say(grid_summary(grid, o));
            if (o.mode == Mode::OrbitalCubes)
                write_orbital_cubes(mol, grid, o, console);
            else
                write_density_cubes(mol, grid, select_density(o, mol, console), o, console);
            break;
        }
        }
        console.finish();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "kumo: " << e.what() << "\nTry 'kumo --help'.\n";
        return 1;
    }
}
