#include "chem/kumo/basis_evaluator.hpp"

#include "basis_math.hpp"
#include "checked_math.hpp"
#include "mo_coefficients.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace chem::kumo {
namespace {

double factorial(int n) {
    if (n < 0) throw std::runtime_error("negative factorial argument");
    double value = 1.0;
    // A wider counter cannot overflow at i == n, even for n == INT_MAX.
    for (long long i = 2; i <= n; ++i)
        value *= static_cast<double>(i);
    return value;
}

double integer_power(double x, int n) noexcept {
    double result = 1.0;
    for (int i = 0; i < n; ++i)
        result *= x;
    return result;
}

// Coefficient of the Cartesian Gaussian x^lx y^ly z^lz (lz = l - lx - ly) in
// the real pure (spherical harmonic) Gaussian with angular momentum l and
// order m, for -l <= m <= l. m = 0 is the zonal function; m > 0 and m < 0 are
// the cosine-like and sine-like combinations, with the sqrt(2) factor of real
// spherical harmonics. Both function types are normalized.
//
// Closed-form expression of H. B. Schlegel and M. J. Frisch, "Transformation
// between Cartesian and pure spherical harmonic Gaussians", Int. J. Quantum
// Chem. 54, 83-87 (1995). In the code:
//   a: ratio of Cartesian to pure normalization constants;
//   b: sqrt((l-|m|)!/(l+|m|)!) / (2^l l!);
//   c: sum over i of binom(l, i) binom(i, j) (-1)^i (2l-2i)! / (l-|m|-2i)!,
//      with j = (l - |m| - lz) / 2, which must be an integer;
//   d, e: sum over k of binom(j, k) binom(|m|, lx-2k) times the sign/sqrt(2)
//      factor that selects the cosine (m > 0) or sine (m < 0) part.
// Returns 0 for components that do not contribute.
double c2pcff(int l, int m, int lx, int ly) {
    const int lz = l - lx - ly;
    const int ma = std::abs(m);
    if (ma > l || lz < 0) return 0.0;

    const double a =
        std::sqrt(factorial(2 * lx) * factorial(2 * ly) * factorial(2 * lz) * factorial(l) /
                  (factorial(2 * l) * factorial(lx) * factorial(ly) * factorial(lz)));
    const double b =
        std::sqrt(factorial(l - ma) / factorial(l + ma)) / (std::pow(2.0, l) * factorial(l));
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
            const double d = factorial(j) / (factorial(k) * factorial(j - k)) * factorial(ma) /
                             (factorial(lx - 2 * k) * factorial(ma + 2 * k - lx));
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

// Radius beyond which every AO of a shell, and each of its first and second
// Cartesian derivatives, is at most `cutoff` in magnitude.
//
// For a primitive with exponent a and Cartesian degree l, each factor of the
// derivative kernel (x^m, m x^(m-1) - 2 a x^(m+1), and the second-derivative
// factor) is bounded by r^(l-2), r^(l-1), ..., r^(l+2) with coefficients
// built from l and a. For r >= 1 every such power is at most r^(l+2), so
//   |AO or derivative| <= transform_norm * sum_p W_p C_p r^(l+2) exp(-a_p r^2),
// with W_p the largest |contraction weight| of primitive p over the shell's
// Cartesian functions, C_p = 1 + l + l^2 + 2 a_p (2l + 2) + 4 a_p^2, and
// transform_norm the largest absolute row sum of the pure transform (1 for
// Cartesian shells). Each term decreases for r > sqrt((l + 2) / (2 a_p)), so
// the bound is found by bisection beyond the largest such radius.
double shell_cutoff_radius(const std::vector<double>& exponents,
                           const Eigen::MatrixXd& primitive_weights, int l, double transform_norm,
                           double cutoff) {
    const auto nprim = exponents.size();
    std::vector<double> coefficient(nprim);
    double start = 1.0;
    for (std::size_t p = 0; p < nprim; ++p) {
        const double a = exponents[p];
        const double weight =
            primitive_weights.col(static_cast<Eigen::Index>(p)).cwiseAbs().maxCoeff();
        coefficient[p] =
            transform_norm * weight * (1.0 + l + l * l + 2.0 * a * (2.0 * l + 2.0) + 4.0 * a * a);
        start = std::max(start, std::sqrt((l + 2.0) / (2.0 * a)));
    }
    const auto bound = [&](double r) {
        double sum = 0.0;
        for (std::size_t p = 0; p < nprim; ++p)
            sum += coefficient[p] * std::pow(r, l + 2) * std::exp(-exponents[p] * r * r);
        return sum;
    };
    if (bound(start) <= cutoff) return start;
    double low = start, high = 2.0 * start;
    while (bound(high) > cutoff)
        high *= 2.0;
    for (int iteration = 0; iteration < 60; ++iteration) {
        const double middle = 0.5 * (low + high);
        (bound(middle) > cutoff ? low : high) = middle;
    }
    return high;
}

} // namespace

// Order of Cartesian components within a shell in Gaussian FCHK files. For
// d and f shells Gaussian uses the explicit orders below (xx, yy, zz, xy, xz,
// yz for d); from g shells on, the order is lexicographic in (lx, ly), i.e.
// zzzz, yzzz, yyzz, yyyz, yyyy, xzzz, ...
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

// (2l+1) x ncart matrix taking a shell's Cartesian components, in Gaussian
// FCHK order, to its pure components in Gaussian order: m = 0, +1, -1, +2,
// -2, ..., +l, -l (row 0 is m = 0, rows 2m-1 and 2m are +m and -m).
Eigen::MatrixXd gaussian_cartesian_to_pure(int l) {
    if (l < 2) throw std::runtime_error("pure transformation requested below d shell");
    const auto gaussian = gaussian_cartesian_order(l);
    Eigen::MatrixXd transform(static_cast<Eigen::Index>(2 * l + 1),
                              static_cast<Eigen::Index>(gaussian.size()));
    transform.setZero();

    for (std::size_t g = 0; g < gaussian.size(); ++g) {
        const auto& p = gaussian[g];
        transform(0, static_cast<Eigen::Index>(g)) = c2pcff(l, 0, p.x, p.y);
        for (int m = 1; m <= l; ++m) {
            transform(2 * m - 1, static_cast<Eigen::Index>(g)) = c2pcff(l, m, p.x, p.y);
            transform(2 * m, static_cast<Eigen::Index>(g)) = c2pcff(l, -m, p.x, p.y);
        }
    }
    return transform;
}

BasisEvaluator::BasisEvaluator(BasisSet basis) : basis_(std::move(basis)) {
    auto normalized = detail::normalized_blocks(basis_);
    blocks_.reserve(normalized.size());
    for (auto& source : normalized) {
        PreparedBlock block;
        block.origin = source.origin;
        block.angular_momentum = source.angular_momentum;
        block.exponents = std::move(source.exponents);
        block.powers = std::move(source.powers);
        block.primitive_weights = std::move(source.weights);
        block.pure_transform = std::move(source.pure_transform);
        block.output_offset = source.offset;
        block.output_count = source.count;
        const double transform_norm =
            block.pure_transform.size() ? block.pure_transform.cwiseAbs().rowwise().sum().maxCoeff()
                                        : 1.0;
        const double radius =
            shell_cutoff_radius(block.exponents, block.primitive_weights, block.angular_momentum,
                                transform_norm, ao_cutoff);
        block.cutoff_radius2 = radius * radius;
        max_cartesian_count_ = std::max(max_cartesian_count_, block.powers.size());
        max_primitive_count_ = std::max(max_primitive_count_, block.exponents.size());
        max_angular_momentum_ = std::max(max_angular_momentum_, block.angular_momentum);
        blocks_.push_back(std::move(block));
    }
}

namespace {
// One-dimensional overlap integral of (x - A)^a exp(-alpha (x - A)^2) and
// (x - B)^b exp(-beta (x - B)^2), evaluated exactly: the Gaussian product
// theorem gives one Gaussian at the product centre, and the polynomial is
// expanded binomially about that centre. Centres in bohr, exponents in
// bohr^-2.
double overlap_axis(int a, int b, double alpha, double beta, double A, double B) {
    const double p = alpha + beta;
    const double center = (alpha * A + beta * B) / p;
    const double PA = center - A, PB = center - B;
    const double prefactor = std::exp(-alpha * beta / p * (A - B) * (A - B));
    std::vector<double> moments(static_cast<std::size_t>(a + b + 1), 0.0);
    moments[0] = std::sqrt(std::numbers::pi / p);
    for (int n = 2; n <= a + b; n += 2)
        moments[static_cast<std::size_t>(n)] =
            (n - 1) * moments[static_cast<std::size_t>(n - 2)] / (2.0 * p);
    double result = 0.0;
    double bin_a = 1.0;
    for (int i = 0; i <= a; ++i) {
        double bin_b = 1.0;
        for (int j = 0; j <= b; ++j) {
            result += bin_a * bin_b * integer_power(PA, a - i) * integer_power(PB, b - j) *
                      moments[static_cast<std::size_t>(i + j)];
            if (j < b) bin_b *= static_cast<double>(b - j) / (j + 1);
        }
        if (i < a) bin_a *= static_cast<double>(a - i) / (i + 1);
    }
    return prefactor * result;
}
} // namespace

Eigen::MatrixXd BasisEvaluator::overlap_matrix() const {
    const Eigen::Index n = static_cast<Eigen::Index>(basis_size());
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
    for (std::size_t a = 0; a < blocks_.size(); ++a) {
        const auto& A = blocks_[a];
        for (std::size_t b = a; b < blocks_.size(); ++b) {
            const auto& B = blocks_[b];
            Eigen::MatrixXd cart =
                Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(A.powers.size()),
                                      static_cast<Eigen::Index>(B.powers.size()));
            for (std::size_t i = 0; i < A.powers.size(); ++i) {
                const auto& x = A.powers[i];
                for (std::size_t j = 0; j < B.powers.size(); ++j) {
                    const auto& y = B.powers[j];
                    double value = 0.0;
                    for (std::size_t p = 0; p < A.exponents.size(); ++p)
                        for (std::size_t q = 0; q < B.exponents.size(); ++q) {
                            const double alpha = A.exponents[p], beta = B.exponents[q];
                            const double primitive =
                                overlap_axis(x.x, y.x, alpha, beta, A.origin.x(), B.origin.x()) *
                                overlap_axis(x.y, y.y, alpha, beta, A.origin.y(), B.origin.y()) *
                                overlap_axis(x.z, y.z, alpha, beta, A.origin.z(), B.origin.z());
                            value += A.primitive_weights(static_cast<Eigen::Index>(i),
                                                         static_cast<Eigen::Index>(p)) *
                                     B.primitive_weights(static_cast<Eigen::Index>(j),
                                                         static_cast<Eigen::Index>(q)) *
                                     primitive;
                        }
                    cart(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = value;
                }
            }
            Eigen::MatrixXd block = cart;
            if (A.pure_transform.size()) block = (A.pure_transform * block).eval();
            if (B.pure_transform.size()) block = (block * B.pure_transform.transpose()).eval();
            const auto ra = static_cast<Eigen::Index>(A.output_offset);
            const auto rb = static_cast<Eigen::Index>(B.output_offset);
            const auto na = static_cast<Eigen::Index>(A.output_count);
            const auto nb = static_cast<Eigen::Index>(B.output_count);
            result.block(ra, rb, na, nb) = block;
            if (a != b) result.block(rb, ra, nb, na) = block.transpose();
        }
    }
    return result;
}

CartesianAOMoments BasisEvaluator::cartesian_moments() const {
    const auto n = static_cast<Eigen::Index>(basis_size());
    CartesianAOMoments result;
    result.overlap = Eigen::MatrixXd::Zero(n, n);
    for (auto& m : result.first)
        m = Eigen::MatrixXd::Zero(n, n);
    for (auto& m : result.second)
        m = Eigen::MatrixXd::Zero(n, n);

    for (std::size_t ia = 0; ia < blocks_.size(); ++ia) {
        const auto& A = blocks_[ia];
        for (std::size_t ib = ia; ib < blocks_.size(); ++ib) {
            const auto& B = blocks_[ib];
            const auto nc = static_cast<Eigen::Index>(A.powers.size());
            const auto mc = static_cast<Eigen::Index>(B.powers.size());
            std::array<Eigen::MatrixXd, 10> cart;
            for (auto& m : cart)
                m = Eigen::MatrixXd::Zero(nc, mc);

            for (Eigen::Index i = 0; i < nc; ++i) {
                const auto& x = A.powers[static_cast<std::size_t>(i)];
                const int ap[3] = {x.x, x.y, x.z};
                for (Eigen::Index j = 0; j < mc; ++j) {
                    const auto& y = B.powers[static_cast<std::size_t>(j)];
                    const int bp[3] = {y.x, y.y, y.z};
                    double sums[10] = {};
                    for (std::size_t p = 0; p < A.exponents.size(); ++p) {
                        for (std::size_t q = 0; q < B.exponents.size(); ++q) {
                            const double alpha = A.exponents[p], beta = B.exponents[q];
                            // Moments about the global origin from overlaps in
                            // which the power of the second (ket) function is
                            // raised: x = (x - B_x) + B_x, so the x moment is
                            // s1 + B_x s0 and the x^2 moment s2 + 2 B_x s1 + B_x^2 s0.
                            double v[3][3];
                            for (int d = 0; d < 3; ++d) {
                                const double c = B.origin[d];
                                const double s0 =
                                    overlap_axis(ap[d], bp[d], alpha, beta, A.origin[d], c);
                                const double s1 =
                                    overlap_axis(ap[d], bp[d] + 1, alpha, beta, A.origin[d], c);
                                const double s2 =
                                    overlap_axis(ap[d], bp[d] + 2, alpha, beta, A.origin[d], c);
                                v[d][0] = s0;
                                v[d][1] = s1 + c * s0;
                                v[d][2] = s2 + 2.0 * c * s1 + c * c * s0;
                            }
                            const double weight =
                                A.primitive_weights(i, static_cast<Eigen::Index>(p)) *
                                B.primitive_weights(j, static_cast<Eigen::Index>(q));
                            const double values[10] = {
                                v[0][0] * v[1][0] * v[2][0], v[0][1] * v[1][0] * v[2][0],
                                v[0][0] * v[1][1] * v[2][0], v[0][0] * v[1][0] * v[2][1],
                                v[0][2] * v[1][0] * v[2][0], v[0][0] * v[1][2] * v[2][0],
                                v[0][0] * v[1][0] * v[2][2], v[0][1] * v[1][1] * v[2][0],
                                v[0][1] * v[1][0] * v[2][1], v[0][0] * v[1][1] * v[2][1]};
                            for (int k = 0; k < 10; ++k)
                                sums[k] += weight * values[k];
                        }
                    }
                    for (int k = 0; k < 10; ++k)
                        cart[k](i, j) = sums[k];
                }
            }
            const auto ra = static_cast<Eigen::Index>(A.output_offset);
            const auto rb = static_cast<Eigen::Index>(B.output_offset);
            const auto na = static_cast<Eigen::Index>(A.output_count);
            const auto nb = static_cast<Eigen::Index>(B.output_count);
            for (int k = 0; k < 10; ++k) {
                Eigen::MatrixXd block = std::move(cart[k]);
                if (A.pure_transform.size()) block = (A.pure_transform * block).eval();
                if (B.pure_transform.size()) block = (block * B.pure_transform.transpose()).eval();
                Eigen::MatrixXd* out =
                    k == 0 ? &result.overlap
                           : (k <= 3 ? &result.first[static_cast<std::size_t>(k - 1)]
                                     : &result.second[static_cast<std::size_t>(k - 4)]);
                out->block(ra, rb, na, nb) = block;
                if (ia != ib) out->block(rb, ra, nb, na) = block.transpose();
            }
        }
    }
    return result;
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
                                      integer_power(r.y(), power.y) * integer_power(r.z(), power.z);
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
        if (r2 > block.cutoff_radius2) {
            output
                .segment(static_cast<Eigen::Index>(block.output_offset),
                         static_cast<Eigen::Index>(block.output_count))
                .setZero();
            continue;
        }
        const auto ncart = block.powers.size();
        const auto nprim = block.exponents.size();

        for (std::size_t p = 0; p < nprim; ++p) {
            radial_scratch[p] = std::exp(-block.exponents[p] * r2);
        }
        // Compute the powers x^n, y^n, z^n once and reuse them for all the
        // Cartesian components of this shell. The primitives are summed in
        // the same order as in evaluate_into_impl_reference().
        const int l = block.angular_momentum;
        double xp_stack[32], yp_stack[32], zp_stack[32];
        std::vector<double> xp_heap, yp_heap, zp_heap;
        double *xp = xp_stack, *yp = yp_stack, *zp = zp_stack;
        if (l >= 32) {
            xp_heap.resize(static_cast<std::size_t>(l + 1));
            yp_heap.resize(static_cast<std::size_t>(l + 1));
            zp_heap.resize(static_cast<std::size_t>(l + 1));
            xp = xp_heap.data();
            yp = yp_heap.data();
            zp = zp_heap.data();
        }
        xp[0] = yp[0] = zp[0] = 1.0;
        for (int n = 1; n <= l; ++n) {
            xp[n] = xp[n - 1] * r.x();
            yp[n] = yp[n - 1] * r.y();
            zp[n] = zp[n - 1] * r.z();
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
    const int panels = derivative_panel_count(order);
    const auto n = static_cast<Eigen::Index>(basis_.function_count());
    Eigen::MatrixXd output(n, panels);
    evaluate_derivatives_batch(point, order, output);
    AODerivatives result;
    result.values = output.col(0);
    result.gradient = output.middleCols(1, 3);
    result.hessian.setZero(n, 6);
    if (order == 2) result.hessian = output.middleCols(4, 6);
    return result;
}

int BasisEvaluator::derivative_panel_count(int order) {
    if (order < 1 || order > 2) throw std::invalid_argument("AO derivative order must be 1 or 2");
    return order == 2 ? 10 : 4;
}

void BasisEvaluator::evaluate_derivatives_batch(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points, int order,
    Eigen::Ref<Eigen::MatrixXd> output) const {
    const int panels = derivative_panel_count(order);
    const Eigen::Index count = points.cols();
    if (output.rows() != static_cast<Eigen::Index>(basis_.function_count()) ||
        output.cols() != panels * count) {
        throw std::runtime_error("AO derivative batch output matrix has wrong dimensions");
    }
    // Scratch reused for every shell and point: primitive exponentials,
    // powers x^m, y^m, z^m, the one-dimensional factors below, and the
    // Cartesian components of pure shells before their transform.
    const auto max_l = static_cast<std::size_t>(max_angular_momentum_);
    const std::size_t degrees = max_l + 1;
    std::vector<double> radial(max_primitive_count_);
    std::vector<double> powers(3 * degrees);
    // factors[((p * 3 + axis) * degrees + m) * 3 + q] holds, for primitive p
    // with exponent a and coordinate x = r[axis], the one-dimensional factor
    // x^m exp(-a x^2) and its first and second derivatives, divided by
    // exp(-a x^2): q = 0: x^m; q = 1: m x^(m-1) - 2 a x^(m+1);
    // q = 2: m (m-1) x^(m-2) - 2 a (2m+1) x^m + 4 a^2 x^(m+2).
    std::vector<double> factors(max_primitive_count_ * 3 * degrees * 3);
    const auto factor = [&](std::size_t p, int axis, int m) {
        return &factors[((p * 3 + static_cast<std::size_t>(axis)) * degrees +
                         static_cast<std::size_t>(m)) *
                        3];
    };
    Eigen::MatrixXd pure_scratch(static_cast<Eigen::Index>(max_cartesian_count_), panels * count);
    Eigen::MatrixXd pure_values(static_cast<Eigen::Index>(max_cartesian_count_), panels * count);
    std::vector<Eigen::Index> inside;
    inside.reserve(static_cast<std::size_t>(count));

    // Shells in the outer loop: each shell fills its rows for all points of
    // the batch, so a pure shell needs a single transform per batch. Only the
    // points within the shell cutoff radius are computed; the others get
    // exact zeros (see ao_cutoff).
    //
    // The primitive Cartesian polynomial and its Gaussian factor are
    // differentiated analytically. The contraction and the pure-shell
    // transform are linear, so they apply unchanged to every derivative.
    for (const auto& block : blocks_) {
        const auto nc = static_cast<Eigen::Index>(block.powers.size());
        const auto nprim = block.exponents.size();
        const int l_shell = block.angular_momentum;
        const bool pure = block.pure_transform.size() != 0;
        const auto offset = static_cast<Eigen::Index>(block.output_offset);
        const auto rows = static_cast<Eigen::Index>(block.output_count);
        auto shell_output = output.middleRows(offset, rows);
        inside.clear();
        for (Eigen::Index j = 0; j < count; ++j)
            if ((points.col(j) - block.origin).squaredNorm() <= block.cutoff_radius2)
                inside.push_back(j);
        const auto n_inside = static_cast<Eigen::Index>(inside.size());
        // Points beyond the cutoff keep these zeros.
        if (n_inside < count) shell_output.setZero();
        if (n_inside == 0) continue;
        // Cartesian shells are written straight to their output columns; pure
        // shells to compact scratch columns (component c of the i-th point
        // within the cutoff at c * n_inside + i), transformed below.
        const auto column = [&](int c, Eigen::Index i) {
            return pure ? c * n_inside + i : c * count + inside[static_cast<std::size_t>(i)];
        };
        using Panel = Eigen::Ref<Eigen::MatrixXd>;
        Panel target =
            pure ? Panel(pure_scratch.topRows(nc)) : Panel(output.middleRows(offset, nc));
        for (Eigen::Index i = 0; i < n_inside; ++i) {
            const Eigen::Vector3d r =
                points.col(inside[static_cast<std::size_t>(i)]) - block.origin;
            const double r2 = r.squaredNorm();
            for (std::size_t p = 0; p < nprim; ++p)
                radial[p] = std::exp(-block.exponents[p] * r2);
            for (int axis = 0; axis < 3; ++axis) {
                double* xp = &powers[static_cast<std::size_t>(axis) * degrees];
                xp[0] = 1.0;
                for (int m = 1; m <= l_shell; ++m)
                    xp[m] = xp[m - 1] * r[axis];
            }
            for (std::size_t p = 0; p < nprim; ++p) {
                const double a = block.exponents[p];
                for (int axis = 0; axis < 3; ++axis) {
                    const double* xp = &powers[static_cast<std::size_t>(axis) * degrees];
                    const double x = r[axis];
                    for (int l = 0; l <= l_shell; ++l) {
                        double* f = factor(p, axis, l);
                        const double xl = xp[l];
                        f[0] = xl;
                        f[1] = (l ? l * xp[l - 1] : 0.0) - 2.0 * a * x * xl;
                        f[2] = order == 2
                                   ? (l >= 2 ? l * (l - 1) * xp[l - 2] : 0.0) -
                                         2.0 * a * (2 * l + 1) * xl + 4.0 * a * a * x * x * xl
                                   : 0.0;
                    }
                }
            }
            for (Eigen::Index k = 0; k < nc; ++k) {
                const auto& power = block.powers[static_cast<std::size_t>(k)];
                double sums[10] = {};
                for (std::size_t p = 0; p < nprim; ++p) {
                    const double weight =
                        block.primitive_weights(k, static_cast<Eigen::Index>(p)) * radial[p];
                    const double* fx = factor(p, 0, power.x);
                    const double* fy = factor(p, 1, power.y);
                    const double* fz = factor(p, 2, power.z);
                    sums[0] += weight * fx[0] * fy[0] * fz[0];
                    sums[1] += weight * fx[1] * fy[0] * fz[0];
                    sums[2] += weight * fx[0] * fy[1] * fz[0];
                    sums[3] += weight * fx[0] * fy[0] * fz[1];
                    if (order == 2) {
                        sums[4] += weight * fx[2] * fy[0] * fz[0];
                        sums[5] += weight * fx[0] * fy[2] * fz[0];
                        sums[6] += weight * fx[0] * fy[0] * fz[2];
                        sums[7] += weight * fx[1] * fy[1] * fz[0];
                        sums[8] += weight * fx[1] * fy[0] * fz[1];
                        sums[9] += weight * fx[0] * fy[1] * fz[1];
                    }
                }
                for (int c = 0; c < panels; ++c)
                    target(k, column(c, i)) = sums[c];
            }
        }
        if (pure) {
            auto transformed = pure_values.topLeftCorner(rows, panels * n_inside);
            transformed.noalias() =
                block.pure_transform * pure_scratch.topLeftCorner(nc, panels * n_inside);
            for (int c = 0; c < panels; ++c)
                for (Eigen::Index i = 0; i < n_inside; ++i)
                    shell_output.col(c * count + inside[static_cast<std::size_t>(i)]) =
                        transformed.col(c * n_inside + i);
        }
    }
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

MolecularOrbitals::MolecularOrbitals(const FormattedCheckpoint& fchk, SpinChannel spin) {
    const auto nbasis = detail::checked_positive_size(fchk.integer("Number of basis functions"),
                                                      "Number of basis functions");
    const char* label =
        spin == SpinChannel::Alpha ? "Alpha MO coefficients" : "Beta MO coefficients";
    if (!fchk.contains(label)) {
        throw std::runtime_error(std::string("missing fchk record: ") + label);
    }
    coefficients_ = detail::read_mo_coefficients(fchk, label, nbasis);
}

std::size_t MolecularOrbitals::basis_size() const noexcept {
    return static_cast<std::size_t>(coefficients_.rows());
}

std::size_t MolecularOrbitals::orbital_count() const noexcept {
    return static_cast<std::size_t>(coefficients_.cols());
}

double MolecularOrbitals::evaluate(std::size_t orbital, const Eigen::VectorXd& ao_values) const {
    if (orbital >= orbital_count()) throw std::out_of_range("MO index out of range");
    if (ao_values.size() != coefficients_.rows())
        throw std::runtime_error("AO vector/MO coefficient size mismatch");
    return coefficients_.col(static_cast<Eigen::Index>(orbital)).dot(ao_values);
}

Eigen::VectorXd MolecularOrbitals::evaluate_all(const Eigen::VectorXd& ao_values) const {
    if (ao_values.size() != coefficients_.rows())
        throw std::runtime_error("AO vector/MO coefficient size mismatch");
    return coefficients_.transpose() * ao_values;
}

} // namespace chem::kumo
