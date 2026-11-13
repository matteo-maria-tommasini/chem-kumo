#pragma once
#include "chem/kumo/density.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace chem::kumo {

/// Hartree in kcal/mol (CODATA 2018), for reporting potentials.
inline constexpr double hartree_in_kcal_per_mol = 627.5094740631;
/// Hartree in eV (CODATA 2018), for reporting energies.
inline constexpr double hartree_in_electronvolt = 27.211386245988;

/// Closed, oriented triangle mesh. Vertex positions are in bohr; each
/// triangle lists its vertices counter-clockwise seen from outside, so that
/// (v1 - v0) x (v2 - v0) points out of the enclosed region.
struct TriangleMesh {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;

    /// Total area, bohr^2.
    [[nodiscard]] double area() const;
    /// Enclosed volume, bohr^3, from the divergence theorem; positive for an
    /// outward-oriented closed mesh.
    [[nodiscard]] double volume() const;
};

/// Statistics of a function V on a surface (Politzer and Murray), from the
/// piecewise-linear interpolation of its vertex values over the triangles:
/// integrals are exact for that interpolation, and the positive and negative
/// regions are bounded by its zero line, which may cross triangles. Values
/// are in the unit of V (hartree/e for the electrostatic potential).
struct SurfaceStatistics {
    double area = 0.0;          ///< A, bohr^2
    double positive_area = 0.0; ///< A+, where V >= 0, bohr^2
    double negative_area = 0.0; ///< A-, where V < 0, bohr^2
    /// Area where |V| <= polar_threshold, bohr^2.
    double nonpolar_area = 0.0;

    double minimum = 0.0;           ///< smallest vertex value
    double maximum = 0.0;           ///< largest vertex value
    std::size_t minimum_vertex = 0; ///< vertex of the minimum
    std::size_t maximum_vertex = 0; ///< vertex of the maximum
    double average = 0.0;           ///< (1/A) integral V
    double positive_average = 0.0;  ///< over A+; NaN if A+ = 0
    double negative_average = 0.0;  ///< over A-; NaN if A- = 0
    double positive_variance = 0.0; ///< sigma^2_+ over A+; 0 if A+ = 0
    double negative_variance = 0.0; ///< sigma^2_- over A-; 0 if A- = 0
    double total_variance = 0.0;    ///< sigma^2_tot = sigma^2_+ + sigma^2_-
    /// Balance of charges nu = sigma^2_+ sigma^2_- / (sigma^2_tot)^2, in
    /// [0, 1/4]; NaN if sigma^2_tot = 0.
    double balance = 0.0;
    /// Internal charge separation Pi = (1/A) integral |V - average|.
    double charge_separation = 0.0;
    /// Molecular polarity index MPI = (1/A) integral |V|.
    double polarity_index = 0.0;

    /// Surface local minima and maxima: vertices whose value is below
    /// (above) that of every neighbour and every neighbour of a neighbour.
    /// Minima are sorted by increasing value, maxima by decreasing value.
    std::vector<std::size_t> local_minima;
    std::vector<std::size_t> local_maxima;
};

/// Statistics of the vertex values `values` (one per mesh vertex) on `mesh`.
/// `polar_threshold` (default 10 kcal/mol in hartree) bounds the nonpolar
/// area. Throws std::invalid_argument if the sizes differ, a value or the
/// threshold is not finite, the threshold is negative, a triangle refers to
/// a missing vertex, or the area is zero.
[[nodiscard]] SurfaceStatistics
surface_statistics(const TriangleMesh& mesh, std::span<const double> values,
                   double polar_threshold = 10.0 / hartree_in_kcal_per_mol);

/// Smallest density isovalue accepted, e/bohr^3. The AO cutoff
/// (BasisEvaluator::ao_cutoff) changes the density by up to about 1e-13
/// e/bohr^3, which would move lower isosurfaces noticeably.
inline constexpr double minimum_density_isovalue = 1e-8;

/// Isosurface rho(r) = isovalue (e/bohr^3) of the total electron density.
///
/// The density is evaluated on an axis-aligned grid of the given spacing
/// (bohr) around the atoms, whose margin grows from 4 bohr until the grid
/// boundary lies outside the isosurface. The grid is split into
/// tetrahedra (six per cube, around its main diagonal), and the surface
/// crosses each edge whose ends lie on either side of the isovalue. Each such
/// vertex is then placed on the analytic isosurface itself, to 1e-10 bohr
/// along its edge (Newton steps safeguarded by bisection). Throws
/// std::invalid_argument for an isovalue below minimum_density_isovalue or
/// not finite, or a nonpositive or nonfinite spacing; std::runtime_error if
/// the isosurface extends 16 bohr beyond the atoms.
[[nodiscard]] TriangleMesh density_isosurface(const Molecule& molecule, double isovalue = 1e-3,
                                              double spacing = 0.25);

/// Options of electrostatic_potential_surface().
struct SurfaceOptions {
    double isovalue = 1e-3; ///< density isosurface, e/bohr^3
    double spacing = 0.25;  ///< grid step of the surface construction, bohr
    /// Error bound of the potential at each vertex, hartree/e (as
    /// EvaluationOptions::potential_tolerance); 0 is exact.
    double potential_tolerance = 1e-8;
};

/// Electrostatic potential on a density isosurface and its statistics.
struct PotentialSurface {
    TriangleMesh mesh;
    std::vector<double> potential; ///< per vertex, hartree/e
    SurfaceStatistics statistics;  ///< polar threshold 10 kcal/mol
};

/// Molecular electrostatic potential (as density::total.electrostatic_potential())
/// on the total-density isosurface of density_isosurface(), with its surface
/// statistics. The result does not depend on the number of OpenMP threads.
/// Throws std::invalid_argument for invalid options.
[[nodiscard]] PotentialSurface electrostatic_potential_surface(const Molecule& molecule,
                                                               const SurfaceOptions& options = {});

/// Average local ionization energy (Sjoberg, Murray, Brinck and Politzer),
/// hartree, at the columns of `points` (bohr):
///
///   I(r) = sum_i n_i |e_i| |phi_i(r)|^2 / sum_i n_i |phi_i(r)|^2,
///
/// over the occupied SCF orbitals phi_i, with energies e_i and occupations
/// n_i (2 for the doubly occupied orbitals of a restricted calculation,
/// otherwise 1), as in Multiwfn. The denominator is the density of the
/// occupied orbitals, which equals the FCHK total density to within the
/// precision of the FCHK matrices (ElectronicStructureInfo). The result does
/// not depend on the number of OpenMP threads.
///
/// Throws std::invalid_argument if a point is not finite, the FCHK has no
/// orbital energies, its occupied MOs do not reproduce its density matrices
/// (ElectronicStructureInfo::ao_mo_eligible(); for example a post-SCF
/// density), or the density vanishes at a point.
[[nodiscard]] std::vector<double> average_local_ionization_energy(
    const Molecule& molecule,
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points);

/// Average local ionization energy on a density isosurface and its statistics.
struct IonizationEnergySurface {
    TriangleMesh mesh;
    std::vector<double> ionization_energy; ///< per vertex, hartree
    /// Statistics of the ionization energy. It is positive, so A+ = A, the
    /// positive average is the average and sigma^2_+ = sigma^2_tot is the
    /// variance about it (Multiwfn's "Variance"); the negative quantities,
    /// nu, the nonpolar area and MPI carry no information.
    SurfaceStatistics statistics;
};

/// Average local ionization energy (average_local_ionization_energy()) on
/// the total-density isosurface of density_isosurface(), with its surface
/// statistics. Throws as density_isosurface() and
/// average_local_ionization_energy().
[[nodiscard]] IonizationEnergySurface local_ionization_energy_surface(const Molecule& molecule,
                                                                      double isovalue = 1e-3,
                                                                      double spacing = 0.25);

} // namespace chem::kumo
