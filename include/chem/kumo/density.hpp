#pragma once
#include "chem/kumo/basis_evaluator.hpp"
#include <Eigen/Core>
#include <filesystem>
#include <functional>
#include <initializer_list>
#include <map>
#include <span>
#include <vector>

namespace chem::kumo {
namespace detail {
struct AoMoTestAccess;
struct SurfaceAccess;
} // namespace detail
/// SCF density component. Difference means alpha minus beta.
/// For MO evaluation only Alpha and Beta are meaningful channels.
enum class Spin { Total, Alpha, Beta, Difference };

/// Default AO screening threshold for grid densities and their derivatives,
/// with both the density-matrix (AO_P) and occupied-MO (AO_MO) contractions.
/// In each block of grid points, AOs whose absolute value (and, for
/// derivatives, that of their derivatives) is at most this at every point
/// are left out of the contraction. The resulting change is at rounding
/// level (about 1e-13 e/bohr^3 in the density in the tested systems).
/// A threshold of 0 disables screening and gives the exact contraction.
inline constexpr double default_screening_threshold = 1e-12;

/// Largest number of grid points per evaluation block (Evaluator and the
/// blocked density methods of Molecule).
inline constexpr std::size_t max_block_size = 65536;

/// Progress report of a grid evaluation: called after each block of grid
/// points with the number of points done so far and the total. Calls are
/// serialized (never concurrent, possibly from an OpenMP worker thread) and
/// `done` increases monotonically up to `total`. An exception thrown by the
/// callback aborts the evaluation and is rethrown to the caller.
using ProgressCallback = std::function<void(std::size_t done, std::size_t total)>;

/// Atomic identity, explicit center charge (in units of e), and Cartesian position (bohr).
/// The center charge is not an atomic partial charge and need not equal Z.
struct Atom {
    int atomic_number = 0;                              ///< element, 1 to 118
    double nuclear_charge = 0.0;                        ///< e; effective charge with ECPs
    Eigen::Vector3d position = Eigen::Vector3d::Zero(); ///< bohr
};

/// Regular 3D grid: r(i,j,k) = origin + steps * [i,j,k].
/// All lengths are in bohr; columns of steps are the grid-axis increments.
/// Values are stored with k fastest, then j, then i.
///
/// Every evaluation function requires a valid grid and otherwise throws
/// std::invalid_argument: positive counts, finite origin and steps, nonzero
/// step vectors, and axes that are not (nearly) linearly dependent. The axis
/// test does not depend on the spacing, so very fine grids are valid.
struct Grid {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Matrix3d steps = Eigen::Matrix3d::Identity(); ///< columns: axis vectors, bohr
    Eigen::Vector3i counts = Eigen::Vector3i::Ones();    ///< points along each axis
    /// Axis-aligned grid with the given spacing over the bounding box of the
    /// atoms enlarged by the margin on every side (at least 2 points per
    /// axis). Throws std::invalid_argument for no atoms, a negative or
    /// nonfinite margin, or a nonpositive spacing.
    static Grid around(const std::vector<Atom>& atoms, double margin_bohr, double spacing_bohr);
    /// Points along the line through `start` and `end`, from `margin_bohr`
    /// before `start` to `margin_bohr` after `end`, evenly spaced with both
    /// ends included: the fewest points (at least 2) whose spacing does not
    /// exceed `spacing_bohr`. Counts are (points, 1, 1); the first axis runs
    /// from start to end and the other two are unit vectors orthogonal to it.
    /// Throws std::invalid_argument for nonfinite or coincident ends, a
    /// negative or nonfinite margin, or a nonpositive spacing.
    static Grid line(const Eigen::Vector3d& start, const Eigen::Vector3d& end, double margin_bohr,
                     double spacing_bohr);
    /// Rectangle in the plane of `a`, `b` and `c`: the first axis runs along
    /// b - a, the second perpendicular to it towards c, and the third is the
    /// unit normal (u x v), with one point. The rectangle covers the three
    /// points enlarged by `margin_bohr` on every side, with both edges of each
    /// axis included and the fewest points whose spacing does not exceed
    /// `spacing_bohr`. Counts are (points along u, points along v, 1).
    /// Throws std::invalid_argument for nonfinite or (nearly) collinear
    /// points, a negative or nonfinite margin, or a nonpositive spacing.
    static Grid plane(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                      double margin_bohr, double spacing_bohr);
};

/// Scalar samples on a grid, with atom positions for Cube output.
/// Density fields have units e/bohr^3; MO amplitude fields bohr^(-3/2).
class Field {
public:
    /// Throws std::invalid_argument if the grid is invalid or the number of
    /// values differs from its number of points.
    Field(Grid grid, std::vector<Atom> atoms, std::vector<double> values);
    /// Read a sample at zero-based grid indices.
    double operator()(int i, int j, int k) const;

    /// All samples in grid order (k fastest).
    std::span<const double> values() const noexcept { return values_; }

    /// Mutable samples without access to vector size or capacity.
    /// The returned view is invalidated by Field destruction or reassignment.
    [[nodiscard]] std::span<double> mutable_values() noexcept { return values_; }

    [[nodiscard]] const Grid& grid() const noexcept { return grid_; }

    [[nodiscard]] const std::vector<Atom>& atoms() const noexcept { return atoms_; }

    /// Write a Gaussian Cube file in bohr, six significant digits per value,
    /// with `comment` as the second header line. Throws std::runtime_error if
    /// the file cannot be written and std::invalid_argument for a nonfinite
    /// value.
    void write_cube(const std::filesystem::path& path,
                    const std::string& comment = "Electron density") const;
    /// Write signed MO amplitudes with Gaussian-style orbital index metadata.
    void write_orbital_cube(const std::filesystem::path& path, std::size_t orbital_index,
                            const std::string& comment = "Molecular orbital amplitude") const;

private:
    Grid grid_;
    std::vector<Atom> atoms_;
    std::vector<double> values_;
};

/// Four collinear SCF density components on the same grid.
struct SpinDensities {
    Field total;
    Field alpha;
    Field beta;
    Field spin;
};

/// Signed MO amplitude and its Cartesian spatial derivatives at one point.
/// Hessian is symmetric; derivatives are in atomic units at fixed geometry.
struct OrbitalDerivatives {
    double value = 0.0;
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
};

/// Spatial derivatives of an electron-density component at one Cartesian point.
/// Units: rho e/bohr^3, gradient e/bohr^4, Hessian/Laplacian e/bohr^5,
/// positive kinetic-energy density tau hartree/bohr^3 (atomic units).
/// All derivatives are with respect to electronic position at fixed geometry.
struct DensityDerivatives {
    double value = 0.0;
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
    double laplacian = 0.0;
    double tau = 0.0;
};
/// Whether the loaded SCF density can be reproduced by occupied MOs.
enum class OrbitalRepresentationStatus { Valid, Inconsistent, Unavailable };

/// Tolerances of the check of the occupied MOs against the FCHK density
/// matrices (ElectronicStructureInfo): the total and the spin density
/// matrices are each reproduced if the Frobenius norm of the difference is
/// at most orbital_check_absolute_tolerance + orbital_check_relative_tolerance
/// times the norm of the total density matrix (the larger of the FCHK and
/// rebuilt ones, and never less than that of the spin matrix).
inline constexpr double orbital_check_relative_tolerance = 1e-8;
inline constexpr double orbital_check_absolute_tolerance = 1e-10;

/// Immutable diagnostic calculated once while loading the FCHK.
struct ElectronicStructureInfo {
    OrbitalRepresentationStatus orbital_status = OrbitalRepresentationStatus::Unavailable;
    /// Differences (Frobenius norm, absolute and relative to the norm of
    /// the total density matrix) between the FCHK total and spin density
    /// matrices and their reconstruction from the occupied MOs; zero when
    /// the status is Unavailable. The spin matrix is rebuilt as a difference
    /// of two matrices of the size of the total one, so its rounding error
    /// (FCHK files store about nine significant digits) scales with the
    /// total matrix; measuring it against the total matrix lets weakly
    /// spin-polarized states, such as broken-symmetry singlets, pass.
    double total_absolute_error = 0.0;
    double total_relative_error = 0.0;
    double spin_absolute_error = 0.0;
    double spin_relative_error = 0.0;
    std::size_t occupied_alpha = 0; ///< occupied alpha orbitals
    std::size_t occupied_beta = 0;  ///< occupied beta orbitals

    /// Whether the occupied MOs reproduce the FCHK density matrices, so that
    /// the AO_MO backend may be used.
    [[nodiscard]] bool ao_mo_eligible() const noexcept {
        return orbital_status == OrbitalRepresentationStatus::Valid;
    }
};

/// Ground-state SCF density and molecular-orbital evaluator from an FCHK.
/// Public calculations use the same stored AO basis and Gaussian ordering.
class Molecule {
public:
    /// Load geometry, basis, SCF density matrices and MOs from a plain or
    /// gzip-compressed FCHK, and check whether the occupied MOs reproduce the
    /// density matrices (see ElectronicStructureInfo). Throws on missing,
    /// malformed or inconsistent records.
    static Molecule from_fchk(const std::filesystem::path& path);

    /// Atoms in FCHK order, positions in bohr.
    const std::vector<Atom>& atoms() const noexcept { return atoms_; }

    /// Result of the AO-MO check done when the FCHK was loaded.
    [[nodiscard]] const ElectronicStructureInfo& electronic_structure_info() const noexcept {
        return electronic_info_;
    }

    /// Evaluate one density component (default: total) in e/bohr^3.
    Field density(const Grid& grid, Spin spin = Spin::Total) const;
    /// Analytic density, spatial gradient, Hessian, Laplacian and positive
    /// kinetic-energy density at one point; all spin components supported.
    [[nodiscard]] DensityDerivatives density_derivatives(const Eigen::Vector3d& point,
                                                         Spin spin = Spin::Total) const;
    /// Evaluate spatial density derivatives at arbitrary Cartesian points.
    /// Each input column is a point in bohr; output order matches the columns.
    /// AO derivatives are shared within each block, with Eigen matrix contractions.
    /// AOs whose value and derivatives are all at most screening_threshold
    /// in a block are left out (see default_screening_threshold); it must be
    /// finite and >= 0, and 0 gives the exact contraction.
    [[nodiscard]] std::vector<DensityDerivatives> density_derivatives_blocked(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        Spin spin = Spin::Total, std::size_t block_size = 128,
        double screening_threshold = default_screening_threshold) const;
    /// Point-by-point reference implementation, kept for numerical regression tests.
    [[nodiscard]] Field density_pointwise(const Grid& grid, Spin spin = Spin::Total) const;
    /// Evaluate total, alpha, beta and alpha-minus-beta together.
    [[nodiscard]] SpinDensities spin_densities(const Grid& grid) const;
    /// Evaluate only requested density components, sharing AO work.
    [[nodiscard]] std::map<Spin, Field> densities(const Grid& grid,
                                                  std::initializer_list<Spin> components) const;
    /// Explicit blocked Eigen path; block_size is in grid points (1 to
    /// max_block_size).
    /// screening_threshold must be finite and >= 0; 0 disables AO screening
    /// (see default_screening_threshold).
    /// `progress`, if set, is called after each block (see ProgressCallback).
    [[nodiscard]] std::map<Spin, Field>
    densities_blocked(const Grid& grid, std::initializer_list<Spin> components,
                      std::size_t block_size = 128,
                      double screening_threshold = default_screening_threshold,
                      const ProgressCallback& progress = {}) const;

    /// Number of AO basis functions.
    std::size_t basis_size() const noexcept { return evaluator_.basis_size(); }

    /// AO coefficient matrix of the FCHK one-particle density operator.
    /// Total and Difference return the stored SCF matrices; Alpha/Beta are
    /// derived from them. No occupied-MO reconstruction is performed.
    /// In a nonorthogonal AO basis, electron count is Tr(P*S).
    [[nodiscard]] Eigen::MatrixXd density_matrix(Spin component = Spin::Total) const;

    /// Whether the FCHK has separate beta MOs (unrestricted SCF).
    [[nodiscard]] bool unrestricted() const noexcept { return unrestricted_; }

    /// Number of available MOs in the selected alpha/beta channel.
    [[nodiscard]] std::size_t orbital_count(Spin channel = Spin::Alpha) const;
    /// Number of occupied spin orbitals (not electron-pair count).
    [[nodiscard]] int occupied_count(Spin channel = Spin::Alpha) const;
    /// MO energies in hartree, ordered consistently with MO indices.
    [[nodiscard]] const std::vector<double>& orbital_energies(Spin channel = Spin::Alpha) const;
    /// Signed MO amplitude; orbital index is 1-based, channel Alpha/Beta.
    [[nodiscard]] Field orbital(const Grid& grid, std::size_t index_one_based,
                                Spin channel = Spin::Alpha) const;
    /// Analytic value, gradient and full Hessian of a signed MO amplitude.
    /// The MO index is 1-based; only Alpha/Beta channels are accepted.
    [[nodiscard]] OrbitalDerivatives orbital_derivatives(const Eigen::Vector3d& point,
                                                         std::size_t index_one_based,
                                                         Spin channel = Spin::Alpha) const;
    /// Evaluate value, gradient and Hessian of one MO at arbitrary points.
    /// Columns of points are Cartesian coordinates in bohr; output order is preserved.
    /// Independent blocks can be evaluated concurrently when OpenMP is enabled.
    [[nodiscard]] std::vector<OrbitalDerivatives> orbital_derivatives_blocked(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        std::size_t index_one_based, Spin channel = Spin::Alpha,
        std::size_t block_size = 128) const;
    /// Evaluate several 1-based MOs with shared AO batches. `progress`, if
    /// set, is called after each block of grid points (see ProgressCallback).
    [[nodiscard]] std::vector<Field> orbitals(const Grid& grid,
                                              const std::vector<std::size_t>& indices_one_based,
                                              Spin channel = Spin::Alpha,
                                              const ProgressCallback& progress = {}) const;

private:
    friend class Evaluator;
    friend struct detail::AoMoTestAccess;
    friend struct detail::SurfaceAccess;
    /// Internal multi-component AO-P contraction; one AO derivative evaluation per point.
    [[nodiscard]] std::map<Spin, std::vector<DensityDerivatives>> density_derivatives_shared(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        const std::vector<Spin>& components, std::size_t block_size, int ao_order = 2,
        double screening_threshold = default_screening_threshold) const;
    [[nodiscard]] std::map<Spin, std::vector<DensityDerivatives>>
    density_derivatives_mo(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                           const std::vector<Spin>& components, std::size_t block_size,
                           int ao_order = 2,
                           double screening_threshold = default_screening_threshold) const;
    Molecule(std::vector<Atom> atoms, BasisEvaluator evaluator, Eigen::MatrixXd total,
             Eigen::MatrixXd spin);
    std::vector<Atom> atoms_;
    BasisEvaluator evaluator_;
    Eigen::MatrixXd total_, spin_;
    Eigen::MatrixXd alpha_mos_, beta_mos_;
    /// MO coefficients of a channel (Alpha or Beta, else throws).
    [[nodiscard]] const Eigen::MatrixXd& mo_coefficients(Spin channel) const;
    std::vector<double> alpha_energies_, beta_energies_;
    int nalpha_ = 0, nbeta_ = 0;
    bool unrestricted_ = false;
    ElectronicStructureInfo electronic_info_;
};
} // namespace chem::kumo
