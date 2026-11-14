#include "cli/summary.hpp"

#include "cli/console.hpp"

#include <iomanip>
#include <sstream>

namespace kumo_cli {
using chem::kumo::Molecule;
using chem::kumo::Spin;

std::string molecule_summary(const Molecule& mol) {
    const auto na = mol.occupied_count(Spin::Alpha), nb = mol.occupied_count(Spin::Beta);
    std::string reference;
    if (mol.unrestricted())
        reference = "unrestricted SCF, " + std::to_string(na) + " alpha and " + std::to_string(nb) +
                    " beta occupied orbitals";
    else if (na == nb)
        reference =
            "restricted closed-shell SCF, " + std::to_string(na) + " doubly occupied orbitals";
    else
        reference = "restricted open-shell SCF, " + std::to_string(na) + " alpha and " +
                    std::to_string(nb) + " beta occupied orbitals";
    return std::to_string(mol.atoms().size()) + " atoms, " + std::to_string(mol.basis_size()) +
           " basis functions; " + reference;
}

// Result of the load-time check of the occupied orbitals against the FCHK
// density matrices, which decides whether AO-MO can be used.
std::string orbital_check_summary(const Molecule& mol) {
    const auto& info = mol.electronic_structure_info();
    if (info.orbital_status == chem::kumo::OrbitalRepresentationStatus::Unavailable)
        return "Occupied orbitals: not available in the FCHK, so AO-MO cannot be used";
    std::ostringstream out;
    out << std::setprecision(2)
        << "Occupied orbitals vs FCHK density matrices: error relative to the total matrix "
        << info.total_relative_error << " (total), " << info.spin_relative_error
        << " (spin); AO-MO needs both <= " << chem::kumo::orbital_check_relative_tolerance;
    if (!info.ao_mo_eligible()) out << " (not met)";
    return out.str();
}

// "Grid: 170 x 178 x 51 = 1,543,260 points, spacing 0.200 bohr, margin 5.000 bohr".
std::string grid_summary(const chem::kumo::Grid& grid, const Options& o) {
    std::ostringstream out;
    out << "Grid: " << grid.counts[0] << " x " << grid.counts[1] << " x " << grid.counts[2] << " = "
        << Console::grouped(static_cast<std::size_t>(grid.counts[0]) *
                            static_cast<std::size_t>(grid.counts[1]) *
                            static_cast<std::size_t>(grid.counts[2]))
        << " points, spacing " << std::fixed << std::setprecision(3) << o.spacing
        << " bohr, margin " << o.margin << " bohr";
    return out.str();
}
} // namespace kumo_cli
