#pragma once
// One-line summaries for --verbose.
#include "chem/kumo.hpp"
#include "cli/options.hpp"
#include <string>

namespace kumo_cli {
/// "60 atoms, 1410 basis functions; restricted closed-shell SCF, 90 doubly
/// occupied orbitals".
std::string molecule_summary(const chem::kumo::Molecule& mol);

/// Result of the load-time check of the occupied orbitals against the FCHK
/// density matrices, which decides whether AO-MO can be used.
std::string orbital_check_summary(const chem::kumo::Molecule& mol);

/// "Grid: 170 x 178 x 51 = 1,543,260 points, spacing 0.200 bohr, margin 5.000 bohr".
std::string grid_summary(const chem::kumo::Grid& grid, const Options& o);
} // namespace kumo_cli
