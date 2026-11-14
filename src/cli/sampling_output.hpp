#pragma once
// Density properties at points, along lines and on planes (--point, --line,
// --plane), as text tables.
#include "chem/kumo.hpp"
#include "cli/console.hpp"
#include "cli/options.hpp"
#include "cli/properties.hpp"

namespace kumo_cli {
/// --point: one table on standard output, one row per point in the order given.
void print_points(const chem::kumo::Molecule& mol, const Options& o,
                  const DensitySelection& selection, Console& console);

/// --line and --plane: evaluate on the line or plane grid and write
/// PREFIX_line.dat or PREFIX_plane.dat, then print the file name.
void write_profile(const chem::kumo::Molecule& mol, const Options& o,
                   const DensitySelection& selection, Console& console);
} // namespace kumo_cli
