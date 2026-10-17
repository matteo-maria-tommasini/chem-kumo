#pragma once
#include "chem/kumo/basis_evaluator.hpp"
#include <Eigen/Core>
#include <filesystem>
#include <initializer_list>
#include <map>
#include <span>
#include <vector>

namespace chem::kumo {
/// SCF density component. Difference means alpha minus beta.
/// For MO evaluation only Alpha and Beta are meaningful channels.
enum class Spin { Total, Alpha, Beta, Difference };
/// Atomic number and Cartesian nuclear position (bohr).
struct Atom { int atomic_number; Eigen::Vector3d position; };
/// Regular 3D grid: r(i,j,k) = origin + steps * [i,j,k].
/// All lengths are in bohr; columns of steps are the grid-axis increments.
struct Grid {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Matrix3d steps = Eigen::Matrix3d::Identity(); // columns: axis vectors, bohr
    Eigen::Vector3i counts = Eigen::Vector3i::Ones();
    static Grid around(const std::vector<Atom>& atoms, double margin_bohr, double spacing_bohr);
};
/// Scalar samples on a grid, with atom positions for Cube output.
/// Density fields have units e/bohr^3; MO amplitude fields bohr^(-3/2).
class Field {
public:
    Field(Grid grid, std::vector<Atom> atoms, std::vector<double> values);
    /// Read a sample at zero-based grid indices.
    double operator()(int i, int j, int k) const;
    std::span<const double> values() const noexcept { return values_; }
    void write_cube(const std::filesystem::path& path, const std::string& comment = "Electron density") const;
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
/// Ground-state SCF density and molecular-orbital evaluator from an FCHK.
/// Public calculations use the same stored AO basis and Gaussian ordering.
class Molecule {
public:
    static Molecule from_fchk(const std::filesystem::path& path);
    const std::vector<Atom>& atoms() const noexcept { return atoms_; }
    /// Evaluate one density component (default: total) in e/bohr^3.
    Field density(const Grid& grid, Spin spin = Spin::Total) const;
    /// Analytic density, spatial gradient, Hessian, Laplacian and positive
    /// kinetic-energy density at one point; all spin components supported.
    [[nodiscard]] DensityDerivatives density_derivatives(
        const Eigen::Vector3d& point, Spin spin = Spin::Total) const;
    /// Evaluate spatial density derivatives at arbitrary Cartesian points.
    /// Each input column is a point in bohr; output order matches the columns.
    /// AO derivatives are shared within each block, with Eigen matrix contractions.
    [[nodiscard]] std::vector<DensityDerivatives> density_derivatives_blocked(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        Spin spin = Spin::Total, std::size_t block_size = 128) const;
    // Reference implementation retained for numerical regression testing.
    [[nodiscard]] Field density_pointwise(const Grid& grid, Spin spin = Spin::Total) const;
    /// Evaluate total, alpha, beta and alpha-minus-beta together.
    [[nodiscard]] SpinDensities spin_densities(const Grid& grid) const;
    /// Evaluate only requested density components, sharing AO work.
    [[nodiscard]] std::map<Spin, Field> densities(const Grid& grid,
        std::initializer_list<Spin> components) const;
    /// Explicit blocked Eigen path; block_size is in grid points.
    [[nodiscard]] std::map<Spin, Field> densities_blocked(const Grid& grid,
        std::initializer_list<Spin> components, std::size_t block_size = 128) const;
    std::size_t basis_size() const noexcept { return evaluator_.basis_size(); }
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
    [[nodiscard]] OrbitalDerivatives orbital_derivatives(
        const Eigen::Vector3d& point, std::size_t index_one_based,
        Spin channel = Spin::Alpha) const;
    /// Evaluate value, gradient and Hessian of one MO at arbitrary points.
    /// Columns of points are Cartesian coordinates in bohr; output order is preserved.
    /// Independent blocks can be evaluated concurrently when OpenMP is enabled.
    [[nodiscard]] std::vector<OrbitalDerivatives> orbital_derivatives_blocked(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        std::size_t index_one_based, Spin channel = Spin::Alpha,
        std::size_t block_size = 128) const;
    /// Evaluate several 1-based MOs with shared AO batches.
    [[nodiscard]] std::vector<Field> orbitals(const Grid& grid,
                                const std::vector<std::size_t>& indices_one_based,
                                Spin channel = Spin::Alpha) const;
private:
    friend class Evaluator;
    // Internal multi-component AO-P contraction; one AO derivative evaluation per point.
    [[nodiscard]] std::map<Spin, std::vector<DensityDerivatives>> density_derivatives_shared(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
        const std::vector<Spin>& components, std::size_t block_size, int ao_order = 2) const;
    Molecule(std::vector<Atom> atoms, BasisEvaluator evaluator, Eigen::MatrixXd total, Eigen::MatrixXd spin);
    std::vector<Atom> atoms_;
    BasisEvaluator evaluator_;
    Eigen::MatrixXd total_, spin_;
    Eigen::MatrixXd alpha_mos_, beta_mos_;
    std::vector<double> alpha_energies_, beta_energies_;
    int nalpha_ = 0, nbeta_ = 0;
    bool unrestricted_ = false;
};
}
