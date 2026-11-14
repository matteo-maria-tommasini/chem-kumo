#pragma once
// Molecular orbitals on the kumo command line (--mo, --list-mo).
#include "chem/kumo.hpp"
#include "cli/console.hpp"
#include "cli/options.hpp"

namespace kumo_cli {
/// --list-mo: the orbitals, occupations and energies of the channels chosen
/// with --spin (alpha, beta or both; default alpha).
void list_orbitals(const chem::kumo::Molecule& mol, const Options& o);

/// --mo: one Cube file per selected orbital of the channels chosen with
/// --spin (alpha, beta or both; default alpha), on `grid`.
void write_orbital_cubes(const chem::kumo::Molecule& mol, const chem::kumo::Grid& grid,
                         const Options& o, Console& console);
} // namespace kumo_cli
