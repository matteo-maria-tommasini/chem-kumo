#pragma once

#include "chem/kumo/basis.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chem::kumo {

struct CartesianPower {
    int x = 0;
    int y = 0;
    int z = 0;
};

[[nodiscard]] std::vector<CartesianPower> gaussian_cartesian_order(int l);
[[nodiscard]] Eigen::MatrixXd gaussian_cartesian_to_pure(int l);

/// AO values and Cartesian spatial derivatives at a fixed nuclear geometry.
/// Columns of gradient are x,y,z; Hessian components are xx,yy,zz,xy,xz,yz.
/// Rows follow the Gaussian FCHK AO ordering. All coordinates are in bohr.
struct AODerivatives {
    Eigen::VectorXd values;
    Eigen::Matrix<double, Eigen::Dynamic, 3> gradient;
    Eigen::Matrix<double, Eigen::Dynamic, 6> hessian;
};

/// Evaluate AO basis functions at points in FCHK AO ordering.
/// Batch output stores one AO vector per column.
class BasisEvaluator {
public:
    explicit BasisEvaluator(BasisSet basis);

    [[nodiscard]] const BasisSet& basis() const noexcept { return basis_; }
    [[nodiscard]] std::size_t basis_size() const noexcept { return basis_.function_count(); }
    [[nodiscard]] Eigen::VectorXd evaluate(const Eigen::Vector3d& point) const;
    /// Exact spatial derivatives through second order, including mixed terms.
    [[nodiscard]] AODerivatives evaluate_derivatives(const Eigen::Vector3d& point) const;
    void evaluate_into(const Eigen::Vector3d& point, Eigen::Ref<Eigen::VectorXd> output) const;
    /// Original pointwise evaluator retained as a regression oracle.
    void evaluate_into_reference(const Eigen::Vector3d& point, Eigen::Ref<Eigen::VectorXd> output) const;
    /// Input: 3 x N point columns; output: nAO x N AO values.
    void evaluate_batch(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                        Eigen::Ref<Eigen::MatrixXd> output) const;

private:
    struct PreparedBlock {
        Eigen::Vector3d origin = Eigen::Vector3d::Zero();
        std::vector<double> exponents;
        std::vector<CartesianPower> powers;
        Eigen::MatrixXd primitive_weights;
        Eigen::MatrixXd pure_transform;
        std::size_t output_offset = 0;
        std::size_t output_count = 0;
    };

    void evaluate_into_impl_reference(const Eigen::Vector3d& point,
                            Eigen::Ref<Eigen::VectorXd> output,
                            std::vector<double>& cartesian_scratch,
                            std::vector<double>& radial_scratch) const;
    void evaluate_into_impl(const Eigen::Vector3d& point,
                            Eigen::Ref<Eigen::VectorXd> output,
                            std::vector<double>& cartesian_scratch,
                            std::vector<double>& radial_scratch) const;

    BasisSet basis_;
    std::vector<PreparedBlock> blocks_;
    std::size_t max_cartesian_count_ = 0;
    std::size_t max_primitive_count_ = 0;
};

enum class OrbitalSpin {
    Alpha,
    Beta,
};

/// MO coefficient matrix evaluator; `orbital` indices here are zero-based.
class MolecularOrbitals {
public:
    explicit MolecularOrbitals(const FormattedCheckpoint& fchk,
                               OrbitalSpin spin = OrbitalSpin::Alpha);

    [[nodiscard]] std::size_t basis_size() const noexcept;
    [[nodiscard]] std::size_t orbital_count() const noexcept;
    [[nodiscard]] double evaluate(std::size_t orbital,
                                  const Eigen::VectorXd& ao_values) const;
    [[nodiscard]] Eigen::VectorXd evaluate_all(const Eigen::VectorXd& ao_values) const;
    [[nodiscard]] const Eigen::MatrixXd& coefficients() const noexcept { return coefficients_; }

private:
    Eigen::MatrixXd coefficients_;
};

} // namespace chem::kumo
