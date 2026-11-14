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

/// Default threshold for leaving negligible atomic orbitals (AOs, the basis
/// functions) out of the evaluation of densities and their derivatives on
/// blocks of points. It applies to both ways of computing a density: from
/// the AO density matrix (AO_P) and from the occupied MOs (AO_MO). In each
/// block, AOs whose absolute value (and, for derivatives, that of their
/// derivatives) is at most this at every point of the block are left out of
/// the sums. The resulting change is at rounding level (about 1e-13 e/bohr^3
/// in the density in the tested systems). A threshold of 0 disables this
/// screening and keeps every AO.
inline constexpr double default_screening_threshold = 1e-12;

/// Largest number of points per evaluation block accepted by Evaluator and
/// by Molecule::densities_blocked().
inline constexpr std::size_t max_block_size = 65536;

/// Progress report of a grid evaluation: called after each block of grid
/// points with the number of points done so far and the total. Calls are
/// serialized (never concurrent, possibly from an OpenMP worker thread) and
/// `done` increases monotonically up to `total`. An exception thrown by the
/// callback aborts the evaluation and is rethrown to the caller.
using ProgressCallback = std::function<void(std::size_t done, std::size_t total)>;

/// An atom: atomic number, nuclear charge (in units of e) and Cartesian
/// position (bohr). The nuclear charge is the one the electrons see; it is
/// smaller than the atomic number when an effective core potential (ECP)
/// replaces core electrons. It is not an atomic partial charge.
struct Atom {
    int atomic_number = 0;                              ///< element, 1 to 118
    double nuclear_charge = 0.0;                        ///< e; effective charge with ECPs
    Eigen::Vector3d position = Eigen::Vector3d::Zero(); ///< bohr
};

/// Regular 3D grid of points r(i, j, k) = origin + steps * (i, j, k), with
/// 0 <= i < counts[0] and likewise for j and k. All lengths are in bohr; the
/// columns of `steps` are the step vectors along the three grid axes, which
/// need not be orthogonal. Values on the grid are stored with k varying
/// fastest, then j, then i ("Field order").
///
/// Every evaluation function requires a valid grid and otherwise throws
/// std::invalid_argument: positive counts, finite origin and steps, nonzero
/// step vectors, and axes that are not (nearly) linearly dependent. The axis
/// test does not depend on the spacing, so very fine grids are valid.
struct Grid {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Matrix3d steps = Eigen::Matrix3d::Identity(); ///< columns: axis vectors, bohr
    Eigen::Vector3i counts = Eigen::Vector3i::Ones();    ///< points along each axis
    /// Grid along the Cartesian axes with step `spacing_bohr`, starting at
    /// the corner of the bounding box of the atoms enlarged by `margin_bohr`
    /// on every side and covering that box (the last point may lie up to one
    /// step beyond it; at least 2 points per axis). Throws
    /// std::invalid_argument for no atoms, an atom with a nonfinite position
    /// or an atomic number below 1, a negative or nonfinite margin, a
    /// nonpositive or nonfinite spacing, or too many points.
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
    /// Rectangle in the plane of `a`, `b` and `c`: the first axis u runs along
    /// b - a, the second v perpendicular to it towards c, and the third is the
    /// unit normal u x v, with one point. The rectangle covers the three
    /// points enlarged by `margin_bohr` on every side, with both edges of each
    /// axis included and the fewest points whose spacing does not exceed
    /// `spacing_bohr`. Counts are (points along u, points along v, 1).
    /// Throws std::invalid_argument for nonfinite or (nearly) collinear
    /// points, a negative or nonfinite margin, or a nonpositive spacing.
    static Grid plane(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                      double margin_bohr, double spacing_bohr);
};

/// Scalar values on a grid, with the atoms (for Cube output). Densities are
/// in e/bohr^3 and MO amplitudes in bohr^(-3/2).
class Field {
public:
    /// Throws std::invalid_argument if the grid is invalid or the number of
    /// values differs from its number of points.
    Field(Grid grid, std::vector<Atom> atoms, std::vector<double> values);
    /// Value at zero-based grid indices; throws std::out_of_range outside
    /// the grid.
    double operator()(int i, int j, int k) const;

    /// All values in Field order (k fastest, see Grid).
    std::span<const double> values() const noexcept { return values_; }

    /// Writable view of the values; their number cannot be changed. The view
    /// becomes invalid when the Field is destroyed or assigned to.
    [[nodiscard]] std::span<double> mutable_values() noexcept { return values_; }

    [[nodiscard]] const Grid& grid() const noexcept { return grid_; }

    [[nodiscard]] const std::vector<Atom>& atoms() const noexcept { return atoms_; }

    /// Write a Gaussian Cube file in bohr, six significant digits per value,
    /// with "chem-kumo <version>" as the first header line and `comment` as
    /// the second. Throws std::runtime_error if the file cannot be written
    /// and std::invalid_argument for a nonfinite value.
    void write_cube(const std::filesystem::path& path,
                    const std::string& comment = "Electron density") const;
    /// Write the field as a Gaussian orbital Cube file: as write_cube(), but
    /// with a negative atom count followed by a line giving one orbital and
    /// its one-based index `orbital_index` (Gaussian's convention). Throws
    /// std::invalid_argument for an index of 0 or above 9999, and as
    /// write_cube().
    void write_orbital_cube(const std::filesystem::path& path, std::size_t orbital_index,
                            const std::string& comment = "Molecular orbital amplitude") const;

private:
    Grid grid_;
    std::vector<Atom> atoms_;
    std::vector<double> values_;
};

/// The four SCF density components on the same grid; `spin` is alpha minus
/// beta.
struct SpinDensities {
    Field total;
    Field alpha;
    Field beta;
    Field spin;
};

/// Signed MO amplitude and its first and second derivatives with respect to
/// the Cartesian coordinates of the point (nuclei fixed), in atomic units.
/// The Hessian is the full symmetric 3 x 3 matrix.
struct OrbitalDerivatives {
    double value = 0.0;
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
};

/// One electron-density component rho and its derivatives at one point,
/// with respect to the Cartesian coordinates of the point (nuclei fixed).
/// Units: rho e/bohr^3, gradient e/bohr^4, Hessian and Laplacian e/bohr^5.
/// tau is the positive kinetic-energy density, in hartree/bohr^3:
/// (1/2) sum_mu,nu P_mu,nu grad phi_mu . grad phi_nu for the AO density
/// matrix P of the component, which for an SCF density equals
/// (1/2) sum_i n_i |grad psi_i|^2 over the occupied MOs psi_i.
struct DensityDerivatives {
    double value = 0.0;
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
    double laplacian = 0.0;
    double tau = 0.0;
};
/// Result of checking, when an FCHK (Gaussian formatted checkpoint file) is
/// loaded, whether its occupied MOs reproduce its SCF density matrices:
/// Valid (they do), Inconsistent (they do not, for example for a post-SCF
/// density), or Unavailable (the FCHK has no or too few MOs).
enum class OrbitalRepresentationStatus { Valid, Inconsistent, Unavailable };

/// Tolerances of the check of the occupied MOs against the FCHK density
/// matrices (see ElectronicStructureInfo): the total and the spin density
/// matrices are each reproduced if the Frobenius norm of the difference is
/// at most orbital_check_absolute_tolerance + orbital_check_relative_tolerance
/// times the norm of the total density matrix (the larger of the FCHK and
/// rebuilt ones, and never less than that of the spin matrix).
inline constexpr double orbital_check_relative_tolerance = 1e-8;
inline constexpr double orbital_check_absolute_tolerance = 1e-10;

/// Result of the check, made once when the FCHK is loaded, of the occupied
/// MOs against the FCHK density matrices.
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
    /// densities may be computed from them (Backend::AO_MO).
    [[nodiscard]] bool ao_mo_eligible() const noexcept {
        return orbital_status == OrbitalRepresentationStatus::Valid;
    }
};

/// Molecule loaded from an FCHK: atoms, basis set, SCF density matrices,
/// MOs and orbital energies. It evaluates ground-state SCF densities, their
/// derivatives and MO amplitudes. The density functions of this class use
/// the FCHK density matrices; Evaluator can also use the occupied MOs.
class Molecule {
public:
    /// Load geometry, basis, SCF density matrices and MOs from a plain or
    /// gzip-compressed FCHK, and check whether the occupied MOs reproduce the
    /// density matrices (see ElectronicStructureInfo). Throws
    /// std::runtime_error if the file cannot be read or has missing,
    /// malformed or inconsistent records.
    static Molecule from_fchk(const std::filesystem::path& path);

    /// Atoms in FCHK order, positions in bohr.
    const std::vector<Atom>& atoms() const noexcept { return atoms_; }

    /// Result of the check of the occupied MOs against the FCHK density
    /// matrices, made when the FCHK was loaded.
    [[nodiscard]] const ElectronicStructureInfo& electronic_structure_info() const noexcept {
        return electronic_info_;
    }

    /// One density component (default: total) on `grid`, in e/bohr^3, with
    /// the default block size and AO screening. Throws std::invalid_argument
    /// for an invalid grid (see Grid).
    Field density(const Grid& grid, DensityComponent spin = DensityComponent::Total) const;
    /// Density, gradient, Hessian, Laplacian and positive kinetic-energy
    /// density of any component at one point (bohr), computed analytically
    /// without AO screening. Throws std::invalid_argument for a nonfinite
    /// point.
    [[nodiscard]] DensityDerivatives
    density_derivatives(const Eigen::Vector3d& point,
                        DensityComponent spin = DensityComponent::Total) const;
    /// density_derivatives() at each column of `points` (bohr), in column
    /// order, evaluated in blocks of `block_size` points with matrix
    /// products. AOs whose value and derivatives are all at most
    /// screening_threshold in a block are left out (see
    /// default_screening_threshold). Throws std::invalid_argument for a zero
    /// block size, a nonfinite point, or a negative or nonfinite threshold.
    [[nodiscard]] std::vector<DensityDerivatives> density_derivatives_blocked(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        DensityComponent spin = DensityComponent::Total, std::size_t block_size = 128,
        double screening_threshold = default_screening_threshold) const;
    /// density() computed point by point, without blocks or AO screening;
    /// kept as a reference for numerical regression tests.
    [[nodiscard]] Field density_pointwise(const Grid& grid,
                                          DensityComponent spin = DensityComponent::Total) const;
    /// All four density components on `grid`, computed together.
    [[nodiscard]] SpinDensities spin_densities(const Grid& grid) const;
    /// The requested density components on `grid`, sharing the AO
    /// evaluation; as densities_blocked() with the default settings.
    [[nodiscard]] std::map<DensityComponent, Field>
    densities(const Grid& grid, std::initializer_list<DensityComponent> components) const;
    /// The requested density components on `grid`, evaluated in blocks of
    /// `block_size` grid points (1 to max_block_size) with AO screening
    /// threshold `screening_threshold` (finite and >= 0; 0 disables it, see
    /// default_screening_threshold). `progress`, if set, is called after each
    /// block (see ProgressCallback). Throws std::invalid_argument for an
    /// invalid grid or an invalid block size or threshold.
    [[nodiscard]] std::map<DensityComponent, Field>
    densities_blocked(const Grid& grid, std::initializer_list<DensityComponent> components,
                      std::size_t block_size = 128,
                      double screening_threshold = default_screening_threshold,
                      const ProgressCallback& progress = {}) const;

    /// Number of AO basis functions.
    std::size_t basis_size() const noexcept { return evaluator_.basis_size(); }

    /// AO density matrix P of one component, so that rho(r) =
    /// sum_mu,nu P_mu,nu phi_mu(r) phi_nu(r). Total and Difference are the
    /// SCF matrices stored in the FCHK; Alpha = (Total + Difference) / 2 and
    /// Beta = (Total - Difference) / 2. The MOs are not used. The number of
    /// electrons is Tr(P S), with S the AO overlap matrix.
    [[nodiscard]] Eigen::MatrixXd
    density_matrix(DensityComponent component = DensityComponent::Total) const;

    /// Whether the FCHK has separate beta MOs (unrestricted SCF).
    [[nodiscard]] bool unrestricted() const noexcept { return unrestricted_; }

    /// Number of MOs of the channel (Alpha or Beta); can be less than the
    /// number of AOs if Gaussian removed linear dependencies.
    [[nodiscard]] std::size_t orbital_count(SpinChannel channel = SpinChannel::Alpha) const;
    /// Number of occupied orbitals of the channel (Alpha or Beta), that is
    /// the number of electrons of that spin.
    [[nodiscard]] std::size_t occupied_count(SpinChannel channel = SpinChannel::Alpha) const;
    /// MO energies of the channel (Alpha or Beta) in hartree; entry i is the
    /// energy of MO i + 1.
    [[nodiscard]] const std::vector<double>&
    orbital_energies(SpinChannel channel = SpinChannel::Alpha) const;
    /// Signed amplitude of MO `index_one_based` (1 = first MO) of the channel
    /// (Alpha or Beta) on `grid`, in bohr^(-3/2). Throws
    /// std::invalid_argument for an invalid grid or channel and
    /// std::out_of_range for an index outside 1 to orbital_count().
    [[nodiscard]] Field orbital(const Grid& grid, std::size_t index_one_based,
                                SpinChannel channel = SpinChannel::Alpha) const;
    /// Value, gradient and Hessian of the signed amplitude of MO
    /// `index_one_based` (1 = first MO) at one point (bohr), computed
    /// analytically. Throws std::invalid_argument for a nonfinite point or a
    /// channel other than Alpha or Beta, and std::out_of_range for an index
    /// outside 1 to orbital_count().
    [[nodiscard]] OrbitalDerivatives
    orbital_derivatives(const Eigen::Vector3d& point, std::size_t index_one_based,
                        SpinChannel channel = SpinChannel::Alpha) const;
    /// orbital_derivatives() at each column of `points` (bohr), in column
    /// order, evaluated in blocks of `block_size` points (in parallel in
    /// OpenMP builds). Throws as orbital_derivatives(), and
    /// std::invalid_argument for a zero block size.
    [[nodiscard]] std::vector<OrbitalDerivatives> orbital_derivatives_blocked(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        std::size_t index_one_based, SpinChannel channel = SpinChannel::Alpha,
        std::size_t block_size = 128) const;
    /// Several MOs of one channel on `grid` (one-based indices), sharing the
    /// AO evaluation; one Field per index, in the order given. `progress`, if
    /// set, is called after each block of grid points (see ProgressCallback).
    /// Throws as orbital(), and std::runtime_error if the FCHK has no MOs.
    [[nodiscard]] std::vector<Field> orbitals(const Grid& grid,
                                              const std::vector<std::size_t>& indices_one_based,
                                              SpinChannel channel = SpinChannel::Alpha,
                                              const ProgressCallback& progress = {}) const;

private:
    friend class Evaluator;
    friend struct detail::AoMoTestAccess;
    friend struct detail::SurfaceAccess;
    /// Derivatives of several density components at the columns of `points`
    /// from the density matrices (AO_P), with one AO evaluation per point
    /// shared by all components; used by Evaluator and the surface functions.
    [[nodiscard]] std::map<DensityComponent, std::vector<DensityDerivatives>>
    density_derivatives_shared(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        const std::vector<DensityComponent>& components, std::size_t block_size, int ao_order = 2,
        double screening_threshold = default_screening_threshold) const;
    [[nodiscard]] std::map<DensityComponent, std::vector<DensityDerivatives>>
    density_derivatives_mo(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                           const std::vector<DensityComponent>& components, std::size_t block_size,
                           int ao_order = 2,
                           double screening_threshold = default_screening_threshold) const;
    Molecule(std::vector<Atom> atoms, BasisEvaluator evaluator, Eigen::MatrixXd total,
             Eigen::MatrixXd spin);
    std::vector<Atom> atoms_;
    BasisEvaluator evaluator_;
    Eigen::MatrixXd total_, spin_;
    Eigen::MatrixXd alpha_mos_, beta_mos_;
    /// MO coefficients of a channel, one MO per column (Alpha or Beta; else
    /// throws std::invalid_argument).
    [[nodiscard]] const Eigen::MatrixXd& mo_coefficients(SpinChannel channel) const;
    std::vector<double> alpha_energies_, beta_energies_;
    int nalpha_ = 0, nbeta_ = 0;
    bool unrestricted_ = false;
    ElectronicStructureInfo electronic_info_;
};
} // namespace chem::kumo
