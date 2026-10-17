#include "chem/kumo/basis_evaluator.hpp"

#include "checked_math.hpp"
#include "basis_math.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace chem::kumo {
namespace {

double factorial(int n) {
    if (n < 0) throw std::runtime_error("negative factorial argument");
    double value = 1.0;
    for (int i = 2; i <= n; ++i) value *= static_cast<double>(i);
    return value;
}

double integer_power(double x, int n) noexcept {
    double result = 1.0;
    for (int i = 0; i < n; ++i) result *= x;
    return result;
}

std::vector<CartesianPower> prism_cartesian_order(int l) {
    std::vector<CartesianPower> result;
    for (int lx = 0; lx <= l; ++lx) {
        for (int ly = 0; ly <= l - lx; ++ly) {
            result.push_back(CartesianPower{lx, ly, l - lx - ly});
        }
    }
    return result;
}

double c2pcff(int l, int m, int lx, int ly) {
    const int lz = l - lx - ly;
    const int ma = std::abs(m);
    if (ma > l || lz < 0) return 0.0;

    const double a = std::sqrt(factorial(2 * lx) * factorial(2 * ly) * factorial(2 * lz) * factorial(l) /
                               (factorial(2 * l) * factorial(lx) * factorial(ly) * factorial(lz)));
    const double b = std::sqrt(factorial(l - ma) / factorial(l + ma)) /
                     (std::pow(2.0, l) * factorial(l));
    const int ilim = (l - ma) / 2;
    const int j = (l - ma - lz) / 2;
    double g = 0.0;
    if (2 * j != l - ma - lz) return 0.0;

    for (int i = 0; i <= ilim; ++i) {
        if (j < 0 || j > i) continue;
        const double c = factorial(l) / (factorial(i) * factorial(l - i)) *
                         factorial(2 * l - 2 * i) * std::pow(-1.0, i) / factorial(l - ma - 2 * i) *
                         factorial(i) / (factorial(j) * factorial(i - j));
        for (int k = 0; k <= j; ++k) {
            if (lx - 2 * k < 0 || lx - 2 * k > ma) continue;
            const double d = factorial(j) / (factorial(k) * factorial(j - k)) *
                             factorial(ma) / (factorial(lx - 2 * k) * factorial(ma + 2 * k - lx));
            double e = 0.0;
            const int k_minus_lx2 = k - lx / 2;
            if (m == 0 && (lx % 2) == 0) e = std::pow(-1.0, k_minus_lx2);
            const int kt = (2 * k + ma - lx) / 2;
            if (m > 0 && (std::abs(ma - lx) % 2) == 0) e = std::sqrt(2.0) * std::pow(-1.0, kt);
            if (m < 0 && (std::abs(ma - lx) % 2) == 1) e = std::sqrt(2.0) * std::pow(-1.0, kt);
            g += c * d * e;
        }
    }
    return a * b * g;
}

} // namespace

std::vector<CartesianPower> gaussian_cartesian_order(int l) {
    if (l < 0) throw std::runtime_error("negative angular momentum");
    if (l == 0) return {{0, 0, 0}};
    if (l == 1) return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    if (l == 2) return {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
    if (l == 3) {
        return {{3, 0, 0}, {0, 3, 0}, {0, 0, 3}, {1, 2, 0}, {2, 1, 0},
                {2, 0, 1}, {1, 0, 2}, {0, 1, 2}, {0, 2, 1}, {1, 1, 1}};
    }

    std::vector<CartesianPower> result;
    for (int lx = 0; lx <= l; ++lx) {
        for (int ly = 0; ly <= l - lx; ++ly) {
            result.push_back(CartesianPower{lx, ly, l - lx - ly});
        }
    }
    return result;
}

Eigen::MatrixXd gaussian_cartesian_to_pure(int l) {
    if (l < 2) throw std::runtime_error("pure transformation requested below d shell");
    const auto gaussian = gaussian_cartesian_order(l);
    const auto prism = prism_cartesian_order(l);
    Eigen::MatrixXd transform(static_cast<Eigen::Index>(2 * l + 1),
                              static_cast<Eigen::Index>(gaussian.size()));
    transform.setZero();

    for (std::size_t g = 0; g < gaussian.size(); ++g) {
        const auto it = std::find_if(prism.begin(), prism.end(), [&](const CartesianPower& p) {
            return p.x == gaussian[g].x && p.y == gaussian[g].y && p.z == gaussian[g].z;
        });
        if (it == prism.end()) throw std::runtime_error("Gaussian/PRISM Cartesian mapping failed");
        const auto& p = *it;
        transform(0, static_cast<Eigen::Index>(g)) = c2pcff(l, 0, p.x, p.y);
        for (int m = 1; m <= l; ++m) {
            transform(2 * m - 1, static_cast<Eigen::Index>(g)) = c2pcff(l, m, p.x, p.y);
            transform(2 * m, static_cast<Eigen::Index>(g)) = c2pcff(l, -m, p.x, p.y);
        }
    }
    return transform;
}

BasisEvaluator::BasisEvaluator(BasisSet basis) : basis_(std::move(basis)) {
    std::size_t output_offset = 0;

    const auto add_block = [&](const Shell& shell,
                               int angular_momentum,
                               const std::vector<double>& coefficients,
                               ShellRepresentation representation,
                               std::vector<PreparedBlock>& blocks,
                               std::size_t& max_cartesian_count,
                               std::size_t& offset) {
        PreparedBlock block;
        block.origin = shell.origin;
        block.exponents = shell.exponents;
        block.powers = gaussian_cartesian_order(angular_momentum);
        block.output_offset = offset;

        const auto ncart = block.powers.size();
        const auto nprim = block.exponents.size();
        block.primitive_weights.resize(static_cast<Eigen::Index>(ncart),
                                       static_cast<Eigen::Index>(nprim));
        for (std::size_t k = 0; k < ncart; ++k) {
            for (std::size_t p = 0; p < nprim; ++p) {
                block.primitive_weights(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(p)) =
                    coefficients[p] * detail::cartesian_primitive_normalization(block.exponents[p], block.powers[k]);
            }
        }

        if (representation == ShellRepresentation::Pure && angular_momentum >= 2) {
            block.pure_transform = gaussian_cartesian_to_pure(angular_momentum);
            block.output_count = static_cast<std::size_t>(block.pure_transform.rows());
        } else {
            block.output_count = ncart;
        }
        max_cartesian_count = std::max(max_cartesian_count, ncart);
        max_primitive_count_ = std::max(max_primitive_count_, nprim);
        offset += block.output_count;
        blocks.push_back(std::move(block));
    };

    blocks_.reserve(basis_.shells().size() + 4U);
    for (const auto& shell : basis_.shells()) {
        if (shell.is_sp) {
            add_block(shell, 0, shell.coefficients, ShellRepresentation::Cartesian,
                      blocks_, max_cartesian_count_, output_offset);
            add_block(shell, 1, shell.p_coefficients, ShellRepresentation::Cartesian,
                      blocks_, max_cartesian_count_, output_offset);
        } else {
            add_block(shell, shell.angular_momentum, shell.coefficients, shell.representation,
                      blocks_, max_cartesian_count_, output_offset);
        }
    }
    if (output_offset != basis_.function_count()) {
        throw std::logic_error("prepared AO evaluator size does not match basis size");
    }
}

void BasisEvaluator::evaluate_into_impl_reference(const Eigen::Vector3d& point,
                                        Eigen::Ref<Eigen::VectorXd> output,
                                        std::vector<double>& cartesian_scratch,
                                        std::vector<double>& radial_scratch) const {
    if (static_cast<std::size_t>(output.size()) != basis_.function_count()) {
        throw std::runtime_error("AO output vector has wrong basis size");
    }
    if (cartesian_scratch.size() < max_cartesian_count_ ||
        radial_scratch.size() < max_primitive_count_) {
        throw std::logic_error("AO evaluator scratch buffer is too small");
    }

    for (const auto& block : blocks_) {
        const Eigen::Vector3d r = point - block.origin;
        const double r2 = r.squaredNorm();
        const auto ncart = block.powers.size();
        const auto nprim = block.exponents.size();

        for (std::size_t p = 0; p < nprim; ++p) {
            radial_scratch[p] = std::exp(-block.exponents[p] * r2);
        }
        for (std::size_t k = 0; k < ncart; ++k) {
            const auto& power = block.powers[k];
            const double polynomial = integer_power(r.x(), power.x) *
                                      integer_power(r.y(), power.y) *
                                      integer_power(r.z(), power.z);
            double radial_sum = 0.0;
            for (std::size_t p = 0; p < nprim; ++p) {
                radial_sum += block.primitive_weights(static_cast<Eigen::Index>(k),
                                                      static_cast<Eigen::Index>(p)) *
                              radial_scratch[p];
            }
            cartesian_scratch[k] = polynomial * radial_sum;
        }

        const auto out_offset = static_cast<Eigen::Index>(block.output_offset);
        if (block.pure_transform.size() == 0) {
            for (std::size_t k = 0; k < ncart; ++k) {
                output[out_offset + static_cast<Eigen::Index>(k)] = cartesian_scratch[k];
            }
        } else {
            const Eigen::Map<const Eigen::VectorXd> cartesian(cartesian_scratch.data(),
                                                               static_cast<Eigen::Index>(ncart));
            output.segment(out_offset, static_cast<Eigen::Index>(block.output_count)).noalias() =
                block.pure_transform * cartesian;
        }
    }
}


void BasisEvaluator::evaluate_into_impl(const Eigen::Vector3d& point,
                                        Eigen::Ref<Eigen::VectorXd> output,
                                        std::vector<double>& cartesian_scratch,
                                        std::vector<double>& radial_scratch) const {
    if (static_cast<std::size_t>(output.size()) != basis_.function_count()) {
        throw std::runtime_error("AO output vector has wrong basis size");
    }
    if (cartesian_scratch.size() < max_cartesian_count_ ||
        radial_scratch.size() < max_primitive_count_) {
        throw std::logic_error("AO evaluator scratch buffer is too small");
    }

    for (const auto& block : blocks_) {
        const Eigen::Vector3d r = point - block.origin;
        const double r2 = r.squaredNorm();
        const auto ncart = block.powers.size();
        const auto nprim = block.exponents.size();

        for (std::size_t p = 0; p < nprim; ++p) {
            radial_scratch[p] = std::exp(-block.exponents[p] * r2);
        }
        // Reuse Cartesian powers across all functions in this shell.
        // Preserve the original primitive summation order for reproducibility.
        const int l = block.powers.empty() ? 0 :
            block.powers.front().x + block.powers.front().y + block.powers.front().z;
        double xp_stack[32], yp_stack[32], zp_stack[32];
        std::vector<double> xp_heap, yp_heap, zp_heap;
        double *xp = xp_stack, *yp = yp_stack, *zp = zp_stack;
        if (l >= 32) {
            xp_heap.resize(static_cast<std::size_t>(l + 1));
            yp_heap.resize(static_cast<std::size_t>(l + 1));
            zp_heap.resize(static_cast<std::size_t>(l + 1));
            xp = xp_heap.data(); yp = yp_heap.data(); zp = zp_heap.data();
        }
        xp[0] = yp[0] = zp[0] = 1.0;
        for (int n = 1; n <= l; ++n) {
            xp[n] = xp[n-1] * r.x();
            yp[n] = yp[n-1] * r.y();
            zp[n] = zp[n-1] * r.z();
        }
        for (std::size_t k = 0; k < ncart; ++k) {
            const auto& power = block.powers[k];
            const double polynomial = (xp[power.x] * yp[power.y]) * zp[power.z];
            double radial_sum = 0.0;
            for (std::size_t p = 0; p < nprim; ++p) {
                radial_sum += block.primitive_weights(static_cast<Eigen::Index>(k),
                                                      static_cast<Eigen::Index>(p)) *
                              radial_scratch[p];
            }
            cartesian_scratch[k] = polynomial * radial_sum;
        }

        const auto out_offset = static_cast<Eigen::Index>(block.output_offset);
        if (block.pure_transform.size() == 0) {
            for (std::size_t k = 0; k < ncart; ++k) {
                output[out_offset + static_cast<Eigen::Index>(k)] = cartesian_scratch[k];
            }
        } else {
            const Eigen::Map<const Eigen::VectorXd> cartesian(cartesian_scratch.data(),
                                                               static_cast<Eigen::Index>(ncart));
            output.segment(out_offset, static_cast<Eigen::Index>(block.output_count)).noalias() =
                block.pure_transform * cartesian;
        }
    }
}

AODerivatives BasisEvaluator::evaluate_derivatives(const Eigen::Vector3d& point) const {
    return evaluate_derivatives(point, 2);
}

AODerivatives BasisEvaluator::evaluate_derivatives(const Eigen::Vector3d& point, int order) const {
    if (order < 1 || order > 2) throw std::invalid_argument("AO derivative order must be 1 or 2");
    const auto n = static_cast<Eigen::Index>(basis_.function_count());
    AODerivatives result;
    result.values = Eigen::VectorXd::Zero(n);
    result.gradient.setZero(n, 3);
    result.hessian.setZero(n, 6);

    // Differentiate the primitive Cartesian polynomial and its Gaussian factor
    // analytically. The contraction and pure-shell transform are linear and
    // therefore applied unchanged to every derivative component.
    for (const auto& block : blocks_) {
        const Eigen::Vector3d r = point - block.origin;
        const double r2 = r.squaredNorm();
        const auto nc = static_cast<Eigen::Index>(block.powers.size());
        Eigen::MatrixXd cart = Eigen::MatrixXd::Zero(nc, order == 2 ? 10 : 4);
        for (Eigen::Index k = 0; k < nc; ++k) {
            const auto& power = block.powers[static_cast<std::size_t>(k)];
            const int degrees[3] = {power.x, power.y, power.z};
            for (std::size_t p = 0; p < block.exponents.size(); ++p) {
                const double a = block.exponents[p];
                const double weight = block.primitive_weights(k, static_cast<Eigen::Index>(p)) *
                                      std::exp(-a * r2);
                double v[3], d[3], dd[3];
                for (int axis = 0; axis < 3; ++axis) {
                    const int l = degrees[axis];
                    const double x = r[axis];
                    const double xl = integer_power(x, l);
                    v[axis] = xl;
                    d[axis] = (l ? l * integer_power(x, l - 1) : 0.0) - 2.0 * a * x * xl;
                    if (order == 2)
                        dd[axis] = (l >= 2 ? l * (l - 1) * integer_power(x, l - 2) : 0.0)
                                 - 2.0 * a * (2 * l + 1) * xl + 4.0 * a * a * x * x * xl;
                }
                cart(k, 0) += weight * v[0] * v[1] * v[2];
                cart(k, 1) += weight * d[0] * v[1] * v[2];
                cart(k, 2) += weight * v[0] * d[1] * v[2];
                cart(k, 3) += weight * v[0] * v[1] * d[2];
                if (order == 2) {
                    cart(k, 4) += weight * dd[0] * v[1] * v[2];
                    cart(k, 5) += weight * v[0] * dd[1] * v[2];
                    cart(k, 6) += weight * v[0] * v[1] * dd[2];
                    cart(k, 7) += weight * d[0] * d[1] * v[2];
                    cart(k, 8) += weight * d[0] * v[1] * d[2];
                    cart(k, 9) += weight * v[0] * d[1] * d[2];
                }
            }
        }
        const auto offset = static_cast<Eigen::Index>(block.output_offset);
        Eigen::MatrixXd transformed;
        if (block.pure_transform.size() == 0) {
            transformed = std::move(cart);
        } else {
            transformed.noalias() = block.pure_transform * cart;
        }
        const auto count = static_cast<Eigen::Index>(block.output_count);
        result.values.segment(offset, count) = transformed.col(0);
        result.gradient.middleRows(offset, count) = transformed.middleCols(1, 3);
        if (order == 2)
            result.hessian.middleRows(offset, count) = transformed.middleCols(4, 6);
    }
    return result;
}

Eigen::VectorXd BasisEvaluator::evaluate(const Eigen::Vector3d& point) const {
    Eigen::VectorXd output(static_cast<Eigen::Index>(basis_.function_count()));
    evaluate_into(point, output);
    return output;
}

void BasisEvaluator::evaluate_into(const Eigen::Vector3d& point,
                                   Eigen::Ref<Eigen::VectorXd> output) const {
    std::vector<double> cartesian_scratch(max_cartesian_count_);
    std::vector<double> radial_scratch(max_primitive_count_);
    evaluate_into_impl(point, output, cartesian_scratch, radial_scratch);
}

void BasisEvaluator::evaluate_into_reference(const Eigen::Vector3d& point,
                                              Eigen::Ref<Eigen::VectorXd> output) const {
    std::vector<double> cartesian_scratch(max_cartesian_count_);
    std::vector<double> radial_scratch(max_primitive_count_);
    evaluate_into_impl_reference(point, output, cartesian_scratch, radial_scratch);
}

void BasisEvaluator::evaluate_batch(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    Eigen::Ref<Eigen::MatrixXd> output) const {
    if (output.rows() != static_cast<Eigen::Index>(basis_.function_count()) ||
        output.cols() != points.cols()) {
        throw std::runtime_error("AO batch output matrix has wrong dimensions");
    }
    std::vector<double> cartesian_scratch(max_cartesian_count_);
    std::vector<double> radial_scratch(max_primitive_count_);
    for (Eigen::Index column = 0; column < points.cols(); ++column) {
        Eigen::Map<Eigen::VectorXd> out_column(output.col(column).data(), output.rows());
        evaluate_into_impl(points.col(column), out_column, cartesian_scratch, radial_scratch);
    }
}

MolecularOrbitals::MolecularOrbitals(const FormattedCheckpoint& fchk, OrbitalSpin spin) {
    const auto nbasis = detail::checked_positive_size(
        fchk.integer("Number of basis functions"), "Number of basis functions");
    const char* label = spin == OrbitalSpin::Alpha ? "Alpha MO coefficients" : "Beta MO coefficients";
    if (!fchk.contains(label)) {
        throw std::runtime_error(std::string("missing fchk record: ") + label);
    }
    const auto& flat = fchk.reals(label);
    if (flat.size() % nbasis != 0U) throw std::runtime_error(std::string("invalid ") + label + " array size");
    const auto nmo = flat.size() / nbasis;
    coefficients_.resize(static_cast<Eigen::Index>(nbasis), static_cast<Eigen::Index>(nmo));
    for (std::size_t mo = 0; mo < nmo; ++mo) {
        for (std::size_t ao = 0; ao < nbasis; ++ao) {
            coefficients_(static_cast<Eigen::Index>(ao), static_cast<Eigen::Index>(mo)) = flat[mo * nbasis + ao];
        }
    }
}

std::size_t MolecularOrbitals::basis_size() const noexcept {
    return static_cast<std::size_t>(coefficients_.rows());
}

std::size_t MolecularOrbitals::orbital_count() const noexcept {
    return static_cast<std::size_t>(coefficients_.cols());
}

double MolecularOrbitals::evaluate(std::size_t orbital, const Eigen::VectorXd& ao_values) const {
    if (orbital >= orbital_count()) throw std::out_of_range("MO index out of range");
    if (ao_values.size() != coefficients_.rows()) throw std::runtime_error("AO vector/MO coefficient size mismatch");
    return coefficients_.col(static_cast<Eigen::Index>(orbital)).dot(ao_values);
}

Eigen::VectorXd MolecularOrbitals::evaluate_all(const Eigen::VectorXd& ao_values) const {
    if (ao_values.size() != coefficients_.rows()) throw std::runtime_error("AO vector/MO coefficient size mismatch");
    return coefficients_.transpose() * ao_values;
}

} // namespace chem::kumo
