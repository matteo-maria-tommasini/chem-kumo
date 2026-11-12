#pragma once

#include "chem/kumo/basis.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <vector>

namespace chem::kumo {

/// Powers (x, y, z) of a Cartesian Gaussian x^x y^y z^z exp(-a r^2).
struct CartesianPower {
    int x = 0;
    int y = 0;
    int z = 0;
};

/// Cartesian components of angular momentum l in Gaussian FCHK order:
/// explicit orders for d and f (xx, yy, zz, xy, xz, yz for d), then
/// lexicographic in (lx, ly) from g on. Throws for negative l.
[[nodiscard]] std::vector<CartesianPower> gaussian_cartesian_order(int l);
/// (2l + 1) x ncart matrix taking normalized Cartesian components, in the
/// order above, to real pure components in Gaussian order m = 0, +1, -1,
/// +2, -2, ... (Schlegel and Frisch 1995). Throws for l < 2.
[[nodiscard]] Eigen::MatrixXd gaussian_cartesian_to_pure(int l);

/// AO values and Cartesian spatial derivatives at a fixed nuclear geometry.
/// Columns of gradient are x,y,z; Hessian components are xx,yy,zz,xy,xz,yz.
/// Rows follow the Gaussian FCHK AO ordering. All coordinates are in bohr.
struct AODerivatives {
    Eigen::VectorXd values;
    Eigen::Matrix<double, Eigen::Dynamic, 3> gradient;
    Eigen::Matrix<double, Eigen::Dynamic, 6> hessian;
};

/// Raw Cartesian AO moments about the global origin, in atomic units.
/// first: x,y,z; second: xx,yy,zz,xy,xz,yz.
/// Matrices follow the Gaussian FCHK AO ordering, including pure shells.
struct CartesianAOMoments {
    Eigen::MatrixXd overlap;
    std::array<Eigen::MatrixXd, 3> first;
    std::array<Eigen::MatrixXd, 6> second;
};

/// Evaluate AO basis functions at points in FCHK AO ordering.
/// Batch output stores one AO vector per column.
class BasisEvaluator {
public:
    /// Shell cutoff: at a point farther from a shell's center than its cutoff
    /// radius, a rigorous bound on the magnitude of each of the shell's AOs
    /// and of their first and second derivatives is at most this value, and
    /// they are returned as exact zeros instead of being computed.
    /// evaluate_into_reference() does not apply the cutoff.
    static constexpr double ao_cutoff = 1e-14;

    /// Precompute normalized primitive weights, pure transforms and shell
    /// cutoff radii.
    explicit BasisEvaluator(BasisSet basis);

    [[nodiscard]] const BasisSet& basis() const noexcept { return basis_; }

    /// Number of AOs.
    [[nodiscard]] std::size_t basis_size() const noexcept { return basis_.function_count(); }

    /// Analytic all-space AO overlap in Gaussian FCHK ordering.
    [[nodiscard]] Eigen::MatrixXd overlap_matrix() const;
    /// Analytic Cartesian moments through order two, computed in one shell-pair traversal.
    [[nodiscard]] CartesianAOMoments cartesian_moments() const;
    /// AO values at one point (bohr), in FCHK order.
    [[nodiscard]] Eigen::VectorXd evaluate(const Eigen::Vector3d& point) const;
    /// Exact spatial derivatives through second order, including mixed terms.
    [[nodiscard]] AODerivatives evaluate_derivatives(const Eigen::Vector3d& point) const;
    /// Compute only derivatives up to order 1 or 2; omitted Hessians are zero.
    [[nodiscard]] AODerivatives evaluate_derivatives(const Eigen::Vector3d& point, int order) const;
    /// AO values at one point into `output`, which must have basis_size()
    /// entries.
    void evaluate_into(const Eigen::Vector3d& point, Eigen::Ref<Eigen::VectorXd> output) const;
    /// Original pointwise evaluator retained as a regression oracle.
    void evaluate_into_reference(const Eigen::Vector3d& point,
                                 Eigen::Ref<Eigen::VectorXd> output) const;
    /// Input: 3 x N point columns; output: nAO x N AO values.
    void evaluate_batch(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                        Eigen::Ref<Eigen::MatrixXd> output) const;
    /// Number of output panels of evaluate_derivatives_batch(): 4 for order 1
    /// (value, x, y, z), 10 for order 2 (followed by xx, yy, zz, xy, xz, yz).
    [[nodiscard]] static int derivative_panel_count(int order);
    /// AO values and derivatives through `order` (1 or 2) at 3 x N points,
    /// without per-point allocations. `output` is nAO x (panels * N), with
    /// panels = derivative_panel_count(order), stored as consecutive nAO x N
    /// panels: component c at point j is column c * N + j. Components follow
    /// AODerivatives: value, gradient x, y, z, then the Hessian in the order
    /// xx, yy, zz, xy, xz, yz. evaluate_derivatives() uses the same kernel.
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
        /// Squared distance from origin beyond which every function of the
        /// shell, and its first and second derivatives, is at most
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

/// MO coefficient set read by MolecularOrbitals.
enum class OrbitalSpin {
    Alpha,
    Beta,
};

/// MO coefficient matrix evaluator; `orbital` indices here are zero-based.
class MolecularOrbitals {
public:
    /// Read the "Alpha MO coefficients" or "Beta MO coefficients" record;
    /// throws std::runtime_error if it is missing or has the wrong size.
    explicit MolecularOrbitals(const FormattedCheckpoint& fchk,
                               OrbitalSpin spin = OrbitalSpin::Alpha);

    [[nodiscard]] std::size_t basis_size() const noexcept;    ///< rows: AOs
    [[nodiscard]] std::size_t orbital_count() const noexcept; ///< columns: MOs
    /// Amplitude of one MO from AO values at a point; throws
    /// std::out_of_range for an index past the last MO and std::runtime_error
    /// for a size mismatch.
    [[nodiscard]] double evaluate(std::size_t orbital, const Eigen::VectorXd& ao_values) const;
    /// Amplitudes of all MOs from AO values at a point.
    [[nodiscard]] Eigen::VectorXd evaluate_all(const Eigen::VectorXd& ao_values) const;

    /// basis_size() x orbital_count() coefficients, one MO per column.
    [[nodiscard]] const Eigen::MatrixXd& coefficients() const noexcept { return coefficients_; }

private:
    Eigen::MatrixXd coefficients_;
};

} // namespace chem::kumo
