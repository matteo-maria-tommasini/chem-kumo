#pragma once
// Cube files of density properties (the default mode of kumo).
#include "chem/kumo.hpp"
#include "cli/console.hpp"
#include "cli/options.hpp"
#include "cli/properties.hpp"

namespace kumo_cli {
/// Evaluate every selected property of every selected component in one
/// Evaluator request and write one Cube file each: first the densities
/// (PREFIX_<component>.cube), then the other properties
/// (PREFIX_<component>_<property>.cube), where PREFIX is the --output prefix
/// and <component> is density (the total), alpha, beta or spin. The
/// gradient is written as its norm. Each file name is printed on standard
/// output.
/// Each Cube comment names the quantity, its units, the component and how it
/// was computed.
void write_density_cubes(const chem::kumo::Molecule& mol, const chem::kumo::Grid& grid,
                         const DensitySelection& selection, const Options& o, Console& console);
} // namespace kumo_cli
