#pragma once

#include "chem/kumo/basis.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <vector>

namespace chem::kumo {

/// Exponents of one Cartesian Gaussian component x^x * y^y * z^z * exp(-a r^2),
/// where the members x, y and z are the powers of the corresponding coordinates.
struct CartesianPower {
    int x = 0;
    int y = 0;
    int z = 0;
};

/// Cartesian components of a shell of angular momentum l, in the order of
/// Gaussian FCHK files: Gaussian's fixed lists for d and f shells (xx, yy,
/// zz, xy, xz, yz for d), and from g shells on increasing (lx, ly) in
/// lexicographic order (zzzz, yzzz, yyzz, ...). Throws std::runtime_error
/// for negative l.
[[nodiscard]] std::vector<CartesianPower> gaussian_cartesian_order(int l);
/// Matrix of size (2l + 1) x (number of Cartesian components) that maps
/// normalized Cartesian components, in the order of
/// gaussian_cartesian_order(), to real pure (spherical) components in
/// Gaussian order m = 0, +1, -1, +2, -2, ... (Schlegel and Frisch 1995).
/// Throws std::runtime_error for l < 2.
[[nodiscard]] Eigen::MatrixXd gaussian_cartesian_to_pure(int l);

/// Values of all atomic orbitals (AOs) at one point and their first and
/// second derivatives with respect to the point's Cartesian coordinates
/// (nuclei fixed). One row per AO, in FCHK order; the gradient columns are
/// x, y, z and the Hessian columns xx, yy, zz, xy, xz, yz. Lengths in bohr.
struct AODerivatives {
    Eigen::VectorXd values;
    Eigen::Matrix<double, Eigen::Dynamic, 3> gradient;
    Eigen::Matrix<double, Eigen::Dynamic, 6> hessian;
};

/// AO moment matrices about the coordinate origin, in atomic units:
/// overlap <mu|nu>, first moments <mu|r_a|nu> (a = x, y, z) and second
/// moments <mu|r_a r_b|nu> (ab = xx, yy, zz, xy, xz, yz), where mu and nu
/// are AOs and r = (x, y, z). Rows and columns are AOs in FCHK order, pure
/// shells included.
struct CartesianAOMoments {
    Eigen::MatrixXd overlap;
    std::array<Eigen::MatrixXd, 3> first;
    std::array<Eigen::MatrixXd, 6> second;
};

/// Evaluates the AOs of a basis set at points (bohr), in FCHK AO order.
/// Batch functions store the AO values of one point per column.
class BasisEvaluator {
public:
    /// Shell cutoff: at a point farther from a shell's center than its cutoff
    /// radius, a rigorous bound on the magnitude of each of the shell's AOs
    /// and of their first and second derivatives is at most this value, and
    /// they are returned as exact zeros instead of being computed.
    /// evaluate_into_reference() does not apply the cutoff.
    static constexpr double ao_cutoff = 1e-14;

    /// Precompute, for every shell, the normalized primitive weights, the
    /// Cartesian-to-pure transform and the cutoff radius (see ao_cutoff).
    explicit BasisEvaluator(BasisSet basis);

    [[nodiscard]] const BasisSet& basis() const noexcept { return basis_; }

    /// Number of AOs.
    [[nodiscard]] std::size_t basis_size() const noexcept { return basis_.function_count(); }

    /// AO overlap matrix S, integrated analytically over all space.
    [[nodiscard]] Eigen::MatrixXd overlap_matrix() const;
    /// Analytic AO overlap, first and second moments, in one pass over shell pairs.
    [[nodiscard]] CartesianAOMoments cartesian_moments() const;
    /// AO values at one point (bohr).
    [[nodiscard]] Eigen::VectorXd evaluate(const Eigen::Vector3d& point) const;
    /// AO values, gradients and Hessians (mixed terms included) at one point.
    [[nodiscard]] AODerivatives evaluate_derivatives(const Eigen::Vector3d& point) const;
    /// AO values and derivatives up to `order` (1 or 2) at one point; with
    /// order 1 the Hessian is returned as zeros. Throws
    /// std::invalid_argument for another order.
    [[nodiscard]] AODerivatives evaluate_derivatives(const Eigen::Vector3d& point, int order) const;
    /// AO values at one point into `output`, which must have basis_size()
    /// entries (otherwise throws std::runtime_error).
    void evaluate_into(const Eigen::Vector3d& point, Eigen::Ref<Eigen::VectorXd> output) const;
    /// Earlier, simpler implementation of evaluate_into(), kept as a
    /// reference for regression tests; it does not apply the AO cutoff.
    void evaluate_into_reference(const Eigen::Vector3d& point,
                                 Eigen::Ref<Eigen::VectorXd> output) const;
    /// AO values at N points: `points` is 3 x N (one point per column, bohr)
    /// and `output` must be nAO x N (nAO = basis_size(); otherwise throws
    /// std::runtime_error).
    void evaluate_batch(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                        Eigen::Ref<Eigen::MatrixXd> output) const;
    /// Number of output panels of evaluate_derivatives_batch(): 4 for order 1
    /// (value, x, y, z), 10 for order 2 (followed by xx, yy, zz, xy, xz, yz).
    /// Throws std::invalid_argument for another order.
    [[nodiscard]] static int derivative_panel_count(int order);
    /// AO values and derivatives up to `order` (1 or 2) at the N columns of
    /// `points`, without allocating per point. `output` must be
    /// nAO x (panels * N), with panels = derivative_panel_count(order)
    /// (otherwise throws std::runtime_error); it holds consecutive nAO x N
    /// panels, so component c at point j is column c * N + j. The components
    /// are ordered as in AODerivatives: value, gradient x, y, z, then the
    /// Hessian xx, yy, zz, xy, xz, yz. evaluate_derivatives() uses the same
    /// code.
    void evaluate_derivatives_batch(
        const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points, int order,
        Eigen::Ref<Eigen::MatrixXd> output) const;

private:
    struct PreparedBlock {
        Eigen::Vector3d origin = Eigen::Vector3d::Zero();
        int angular_momentum = 0;
        std::vector<double> exponents;
        std::vector<CartesianPower> powers;
        Eigen::MatrixXd primitive_weights;
        Eigen::MatrixXd pure_transform;
        std::size_t output_offset = 0;
        std::size_t output_count = 0;
        /// Squared distance from the shell centre beyond which every AO of
        /// the shell, and its first and second derivatives, is at most
        /// ao_cutoff in magnitude; such points get exact zeros.
        double cutoff_radius2 = 0.0;
    };

    void evaluate_into_impl_reference(const Eigen::Vector3d& point,
                                      Eigen::Ref<Eigen::VectorXd> output,
                                      std::vector<double>& cartesian_scratch,
                                      std::vector<double>& radial_scratch) const;
    void evaluate_into_impl(const Eigen::Vector3d& point, Eigen::Ref<Eigen::VectorXd> output,
                            std::vector<double>& cartesian_scratch,
                            std::vector<double>& radial_scratch) const;

    BasisSet basis_;
    std::vector<PreparedBlock> blocks_;
    std::size_t max_cartesian_count_ = 0;
    std::size_t max_primitive_count_ = 0;
    int max_angular_momentum_ = 0;
};

/// Which molecular-orbital (MO) coefficient record MolecularOrbitals reads.
enum class OrbitalSpin {
    Alpha,
    Beta,
};

/// MO coefficient matrix of one spin, and MO amplitudes computed from AO
/// values. Unlike Molecule, `orbital` indices here are zero-based.
class MolecularOrbitals {
public:
    /// Read the "Alpha MO coefficients" or "Beta MO coefficients" record;
    /// throws std::runtime_error if it is missing or has the wrong size.
    explicit MolecularOrbitals(const FormattedCheckpoint& fchk,
                               OrbitalSpin spin = OrbitalSpin::Alpha);

    [[nodiscard]] std::size_t basis_size() const noexcept;    ///< rows: AOs
    [[nodiscard]] std::size_t orbital_count() const noexcept; ///< columns: MOs
    /// Amplitude of one MO at a point, from the AO values there; throws
    /// std::out_of_range for an index past the last MO and std::runtime_error
    /// for a size mismatch.
    [[nodiscard]] double evaluate(std::size_t orbital, const Eigen::VectorXd& ao_values) const;
    /// Amplitudes of all MOs at a point, from the AO values there; throws
    /// std::runtime_error for a size mismatch.
    [[nodiscard]] Eigen::VectorXd evaluate_all(const Eigen::VectorXd& ao_values) const;

    /// basis_size() x orbital_count() coefficients, one MO per column.
    [[nodiscard]] const Eigen::MatrixXd& coefficients() const noexcept { return coefficients_; }

private:
    Eigen::MatrixXd coefficients_;
};

} // namespace chem::kumo
