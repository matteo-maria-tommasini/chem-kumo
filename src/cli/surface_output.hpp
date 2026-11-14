#pragma once
// Surface reports (--esp-surface, --alie-surface, --surface-map) and PLY
// files (--ply).
#include "chem/kumo.hpp"
#include "cli/console.hpp"
#include "cli/options.hpp"
#include "cli/properties.hpp"

namespace kumo_cli {
/// --esp-surface: statistics of the electrostatic potential on the density
/// isosurface, as a report on standard output.
void report_esp_surface(const chem::kumo::Molecule& mol, const Options& o, Console& console);

/// --alie-surface: statistics of the average local ionization energy on the
/// density isosurface, as a report on standard output.
void report_alie_surface(const chem::kumo::Molecule& mol, const Options& o, Console& console);

/// --surface-map: statistics of each selected property on the density
/// isosurface.
void report_surface_map(const chem::kumo::Molecule& mol, const Options& o,
                        const DensitySelection& selection, Console& console);
} // namespace kumo_cli
