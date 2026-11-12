#include "electrostatic_potential.hpp"
#include "basis_math.hpp"
#include "boys_function.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace chem::kumo::detail {
namespace {
// One block of AOs with a single angular momentum, laid out exactly as in
// BasisEvaluator: SP shells become an s block followed by a p block.
struct Block {
    Eigen::Vector3d origin;
    int angular_momentum = 0;
    std::vector<double> exponents;
    std::vector<CartesianPower> powers;
    Eigen::MatrixXd weights;   ///< Cartesian component x primitive, normalized
    Eigen::MatrixXd transform; ///< pure x Cartesian, or empty for Cartesian
    Eigen::Index offset = 0;
    Eigen::Index count = 0;
};

std::vector<Block> blocks_of(const BasisSet& basis) {
    std::vector<Block> blocks;
    Eigen::Index offset = 0;
    const auto add = [&](const Shell& shell, int l, const std::vector<double>& coefficients,
                         ShellRepresentation representation) {
        if (l > boys_max_order / 2)
            throw std::invalid_argument("the electrostatic potential supports shells up to l = " +
                                        std::to_string(boys_max_order / 2));
        Block block;
        block.origin = shell.origin;
        block.angular_momentum = l;
        block.exponents = shell.exponents;
        block.powers = gaussian_cartesian_order(l);
        const auto ncart = static_cast<Eigen::Index>(block.powers.size());
        const auto nprim = static_cast<Eigen::Index>(block.exponents.size());
        block.weights.resize(ncart, nprim);
        for (Eigen::Index k = 0; k < ncart; ++k)
            for (Eigen::Index p = 0; p < nprim; ++p)
                block.weights(k, p) =
                    coefficients[static_cast<std::size_t>(p)] *
                    cartesian_primitive_normalization(block.exponents[static_cast<std::size_t>(p)],
                                                      block.powers[static_cast<std::size_t>(k)]);
        if (representation == ShellRepresentation::Pure && l >= 2)
            block.transform = gaussian_cartesian_to_pure(l);
        block.offset = offset;
        block.count = block.transform.size() ? block.transform.rows() : ncart;
        offset += block.count;
        blocks.push_back(std::move(block));
    };
    for (const auto& shell : basis.shells()) {
        if (shell.is_sp) {
            add(shell, 0, shell.coefficients, ShellRepresentation::Cartesian);
            add(shell, 1, shell.p_coefficients, ShellRepresentation::Cartesian);
        } else {
            add(shell, shell.angular_momentum, shell.coefficients, shell.representation);
        }
    }
    if (static_cast<std::size_t>(offset) != basis.function_count())
        throw std::logic_error("electrostatic potential blocks do not match the basis size");
    return blocks;
}

// Hermite expansion coefficients E^(ij)_t of the one-dimensional product
// (x - A)^i (x - B)^j exp(-alpha (x-A)^2 - beta (x-B)^2), without the
// exponential prefactor exp(-mu (A - B)^2), which is applied to the whole
// three-dimensional pair:
//   E^(i+1,j)_t = E^(ij)_(t-1) / (2p) + X_PA E^(ij)_t + (t + 1) E^(ij)_(t+1),
//   E^(i,j+1)_t = E^(ij)_(t-1) / (2p) + X_PB E^(ij)_t + (t + 1) E^(ij)_(t+1).
// Stored as table[(i * (lb + 1) + j) * (la + lb + 1) + t].
class HermiteTable {
public:
    HermiteTable(int la, int lb, double p, double pa, double pb)
        : lb_(lb), width_(la + lb + 1),
          values_(static_cast<std::size_t>((la + 1) * (lb + 1) * (la + lb + 1)), 0.0) {
        const double half = 0.5 / p;
        at(0, 0, 0) = 1.0;
        const auto step = [&](int si, int sj, int ti, int tj, double shift) {
            for (int t = 0; t <= ti + tj; ++t) {
                double value = shift * get(si, sj, t);
                if (t > 0) value += half * get(si, sj, t - 1);
                if (t + 1 <= si + sj) value += (t + 1) * get(si, sj, t + 1);
                at(ti, tj, t) = value;
            }
        };
        for (int i = 0; i < la; ++i)
            step(i, 0, i + 1, 0, pa);
        for (int i = 0; i <= la; ++i)
            for (int j = 0; j < lb; ++j)
                step(i, j, i, j + 1, pb);
    }

    [[nodiscard]] double get(int i, int j, int t) const {
        if (t < 0 || t > i + j) return 0.0;
        return values_[index(i, j, t)];
    }

private:
    [[nodiscard]] std::size_t index(int i, int j, int t) const {
        return static_cast<std::size_t>((i * (lb_ + 1) + j) * width_ + t);
    }

    double& at(int i, int j, int t) { return values_[index(i, j, t)]; }

    int lb_, width_;
    std::vector<double> values_;
};

struct HermiteIndex {
    int t, u, v;
};

// (t, u, v) with t + u + v <= L, ordered by t + u + v.
std::vector<HermiteIndex> hermite_triples(int L) {
    std::vector<HermiteIndex> triples;
    for (int sum = 0; sum <= L; ++sum)
        for (int t = sum; t >= 0; --t)
            for (int u = sum - t; u >= 0; --u)
                triples.push_back({t, u, sum - t - u});
    return triples;
}
} // namespace

struct ElectrostaticPotential::Scratch {
    std::vector<double> boys;
    std::vector<double> previous, current;
    std::vector<HermiteIndex> triples;
    std::vector<double> local; ///< local expansion coefficients
};

ElectrostaticPotential::ElectrostaticPotential(const BasisSet& basis,
                                               const Eigen::MatrixXd& density_matrix,
                                               std::vector<Atom> atoms)
    : atoms_(std::move(atoms)) {
    const auto n = static_cast<Eigen::Index>(basis.function_count());
    if (density_matrix.rows() != n || density_matrix.cols() != n)
        throw std::invalid_argument("density matrix does not match the basis");
    if (!density_matrix.allFinite()) throw std::invalid_argument("density matrix must be finite");
    for (const auto& atom : atoms_)
        if (!atom.position.allFinite() || !std::isfinite(atom.nuclear_charge))
            throw std::invalid_argument("nuclear positions and charges must be finite");

    const auto blocks = blocks_of(basis);
    for (const auto& block : blocks)
        max_order_ = std::max(max_order_, 2 * block.angular_momentum);
    side_ = max_order_ + local_order + 1;
    const int side = side_;
    hermite_indices_.resize(static_cast<std::size_t>(max_order_ + 1));
    std::vector<std::vector<HermiteIndex>> triples(static_cast<std::size_t>(max_order_ + 1));
    for (int L = 0; L <= max_order_; ++L) {
        triples[static_cast<std::size_t>(L)] = hermite_triples(L);
        for (const auto& h : triples[static_cast<std::size_t>(L)])
            hermite_indices_[static_cast<std::size_t>(L)].push_back((h.t * side + h.u) * side +
                                                                    h.v);
    }

    // Local expansion terms beta, |beta| <= local_order: their offsets in the
    // Hermite table (added to those of alpha), (-1)^|beta| / beta!, and the
    // powers of h.
    for (const auto& h : hermite_triples(local_order)) {
        local_offsets_.push_back((h.t * side + h.u) * side + h.v);
        double factor = ((h.t + h.u + h.v) % 2 == 0) ? 1.0 : -1.0;
        for (int j = 2; j <= h.t; ++j)
            factor /= j;
        for (int j = 2; j <= h.u; ++j)
            factor /= j;
        for (int j = 2; j <= h.v; ++j)
            factor /= j;
        local_factors_.push_back(factor);
        local_powers_.push_back({h.t, h.u, h.v});
    }

    std::vector<double> pair_coefficients;
    for (std::size_t a = 0; a < blocks.size(); ++a) {
        const Block& A = blocks[a];
        for (std::size_t b = a; b < blocks.size(); ++b) {
            const Block& B = blocks[b];
            // Cartesian density block; the (b, a) block is its transpose and
            // gives the same contribution, hence the factor 2.
            Eigen::MatrixXd cartesian = density_matrix.block(A.offset, B.offset, A.count, B.count);
            if (A.transform.size()) cartesian = (A.transform.transpose() * cartesian).eval();
            if (B.transform.size()) cartesian = (cartesian * B.transform).eval();
            if (a != b) cartesian *= 2.0;
            if ((cartesian.array() == 0.0).all()) continue;

            const int la = A.angular_momentum, lb = B.angular_momentum, L = la + lb;
            const auto& list = triples[static_cast<std::size_t>(L)];
            const Eigen::Vector3d AB = A.origin - B.origin;
            for (std::size_t p = 0; p < A.exponents.size(); ++p) {
                for (std::size_t q = 0; q < B.exponents.size(); ++q) {
                    const double alpha = A.exponents[p], beta = B.exponents[q];
                    const double exponent = alpha + beta;
                    const Eigen::Vector3d center = (alpha * A.origin + beta * B.origin) / exponent;
                    const double scale = std::exp(-alpha * beta / exponent * AB.squaredNorm()) *
                                         2.0 * std::numbers::pi / exponent;
                    const HermiteTable ex(la, lb, exponent, center.x() - A.origin.x(),
                                          center.x() - B.origin.x());
                    const HermiteTable ey(la, lb, exponent, center.y() - A.origin.y(),
                                          center.y() - B.origin.y());
                    const HermiteTable ez(la, lb, exponent, center.z() - A.origin.z(),
                                          center.z() - B.origin.z());
                    pair_coefficients.assign(list.size(), 0.0);
                    for (std::size_t i = 0; i < A.powers.size(); ++i) {
                        const auto& pi = A.powers[i];
                        for (std::size_t j = 0; j < B.powers.size(); ++j) {
                            const auto& pj = B.powers[j];
                            const double weight = cartesian(static_cast<Eigen::Index>(i),
                                                            static_cast<Eigen::Index>(j)) *
                                                  A.weights(static_cast<Eigen::Index>(i),
                                                            static_cast<Eigen::Index>(p)) *
                                                  B.weights(static_cast<Eigen::Index>(j),
                                                            static_cast<Eigen::Index>(q));
                            if (weight == 0.0) continue;
                            for (std::size_t k = 0; k < list.size(); ++k) {
                                const auto& h = list[k];
                                if (h.t > pi.x + pj.x || h.u > pi.y + pj.y || h.v > pi.z + pj.z)
                                    continue;
                                pair_coefficients[k] += weight * ex.get(pi.x, pj.x, h.t) *
                                                        ey.get(pi.y, pj.y, h.u) *
                                                        ez.get(pi.z, pj.z, h.v);
                            }
                        }
                    }
                    if (std::all_of(pair_coefficients.begin(), pair_coefficients.end(),
                                    [](double c) { return c == 0.0; }))
                        continue;
                    pairs_.push_back({exponent, center, L, coefficients_.size()});
                    for (const double c : pair_coefficients)
                        coefficients_.push_back(scale * c);
                }
            }
        }
    }
}

double ElectrostaticPotential::nuclear(const Eigen::Vector3d& point) const {
    if (!point.allFinite())
        throw std::invalid_argument("electrostatic potential point must be finite");
    double sum = 0.0;
    for (const auto& atom : atoms_) {
        const double distance = (atom.position - point).norm();
        // On a nucleus its own term is singular and is left out (see the
        // class documentation).
        if (atom.nuclear_charge == 0.0 || distance < nucleus_exclusion_radius) continue;
        sum += atom.nuclear_charge / distance;
    }
    return sum;
}

double ElectrostaticPotential::electronic(const Eigen::Vector3d& point) const {
    Scratch scratch;
    return electronic(point, scratch);
}

void ElectrostaticPotential::prepare(Scratch& scratch) const {
    const auto cube = static_cast<std::size_t>(side_ * side_ * side_);
    if (scratch.current.size() == cube) return;
    scratch.boys.assign(static_cast<std::size_t>(side_), 0.0);
    scratch.previous.assign(cube, 0.0);
    scratch.current.assign(cube, 0.0);
    scratch.triples = hermite_triples(side_ - 1);
    scratch.local.assign(local_offsets_.size(), 0.0);
}

// Hermite Coulomb integrals (Helgaker, Jorgensen and Olsen, Molecular
// Electronic-Structure Theory, section 9.9):
//   R^n_000 = (-2p)^n F_n(p |PC|^2),
//   R^n_(t+1,u,v) = t R^(n+1)_(t-1,u,v) + X_PC R^(n+1)_(t,u,v),
// and likewise for u and v, from n = L down to n = 0. R^0_tuv for
// t + u + v <= L is left in scratch.current, at offset (t * side_ + u) *
// side_ + v. The integral of a Hermite Gaussian Lambda_tuv over 1/|r - C|
// is (2 pi / p) R^0_tuv; the factor 2 pi / p is in the stored coefficients.
void ElectrostaticPotential::hermite_integrals(double exponent, const Eigen::Vector3d& pc, int L,
                                               Scratch& scratch) const {
    const int side = side_;
    const auto offset = [side](const HermiteIndex& h) {
        return static_cast<std::size_t>((h.t * side + h.u) * side + h.v);
    };
    const std::span<double> boys(scratch.boys.data(), static_cast<std::size_t>(L + 1));
    boys_function(exponent * pc.squaredNorm(), boys);
    // Triples with t + u + v <= s are the first C(s + 3, 3) entries.
    const auto count = [](int s) {
        return static_cast<std::size_t>((s + 1) * (s + 2) * (s + 3) / 6);
    };
    double power = 1.0;
    for (int n = 0; n < L; ++n)
        power *= -2.0 * exponent;
    // Level n holds R^n for t + u + v <= L - n.
    for (int n = L; n >= 0; --n) {
        std::swap(scratch.previous, scratch.current);
        auto& current = scratch.current;
        const auto& previous = scratch.previous;
        current[0] = power * boys[static_cast<std::size_t>(n)];
        const std::size_t entries = count(L - n);
        for (std::size_t k = 1; k < entries; ++k) {
            const auto& h = scratch.triples[k];
            double value;
            if (h.t > 0) {
                value = pc.x() * previous[offset({h.t - 1, h.u, h.v})];
                if (h.t > 1) value += (h.t - 1) * previous[offset({h.t - 2, h.u, h.v})];
            } else if (h.u > 0) {
                value = pc.y() * previous[offset({h.t, h.u - 1, h.v})];
                if (h.u > 1) value += (h.u - 1) * previous[offset({h.t, h.u - 2, h.v})];
            } else {
                value = pc.z() * previous[offset({h.t, h.u, h.v - 1})];
                if (h.v > 1) value += (h.v - 1) * previous[offset({h.t, h.u, h.v - 2})];
            }
            current[offset(h)] = value;
        }
        if (n > 0) power /= -2.0 * exponent;
    }
}

double ElectrostaticPotential::pair_potential(const PrimitivePair& pair,
                                              const Eigen::Vector3d& point,
                                              Scratch& scratch) const {
    hermite_integrals(pair.exponent, pair.center - point, pair.order, scratch);
    const auto& indices = hermite_indices_[static_cast<std::size_t>(pair.order)];
    const double* c = coefficients_.data() + pair.first_coefficient;
    double sum = 0.0;
    for (std::size_t k = 0; k < indices.size(); ++k)
        sum += c[k] * scratch.current[static_cast<std::size_t>(indices[k])];
    return sum;
}

double ElectrostaticPotential::electronic(const Eigen::Vector3d& point, Scratch& scratch) const {
    if (!point.allFinite())
        throw std::invalid_argument("electrostatic potential point must be finite");
    prepare(scratch);
    double sum = 0.0;
    for (const auto& pair : pairs_)
        sum += pair_potential(pair, point, scratch);
    return -sum;
}

void ElectrostaticPotential::evaluate(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    std::span<double> output) const {
    if (output.size() != static_cast<std::size_t>(points.cols()))
        throw std::invalid_argument(
            "electrostatic potential output size does not match the points");
    Scratch scratch;
    for (Eigen::Index j = 0; j < points.cols(); ++j) {
        const Eigen::Vector3d point = points.col(j);
        output[static_cast<std::size_t>(j)] = nuclear(point) + electronic(point, scratch);
    }
}

// Error bound of the local expansion of one primitive pair (see the class
// documentation). Its potential is
//   sum_k c_k R_(alpha_k)(p, P - C) = sum_k e_k d^alpha_k [1 / R]
//                                     - sum_k e_k d^alpha_k [erfc(sqrt(p) R) / R],
// with R = |P - C| and e_k = c_k sqrt(pi / p) / 2. For a cluster of points
// within rho of its centre b, at distance d = |P - b| > rho, the error of the
// Taylor polynomial of order K about b is at most
//   sum_L W_L [ T_L(d) + G_L(d - rho) + sum_(n<=K) (sqrt(3) rho)^n / n! G_(L+n)(d) ],
// with W_L the sum of |e_k| over the terms of order L:
//  - T_L(d) = sum_(n>K) (n+L)! / n! rho^n / d^(n+L+1) bounds the Taylor
//    remainder of the 1/R terms: every derivative of order m of 1/|x| has
//    norm at most m! / |x|^(m+1) as a symmetric multilinear form (Legendre
//    polynomials are bounded by 1, and by Banach's theorem the norm of a
//    symmetric form on a Hilbert space equals that of its diagonal);
//  - G_m(r) bounds every derivative of order m of erfc(sqrt(p) R) / R at
//    R >= r, from erfc(sqrt(p) R) / R = (2/sqrt(pi)) int_sqrt(p)^inf
//    exp(-s^2 R^2) ds and Cramer's inequality |H_n(x)| exp(-x^2/2) <=
//    k sqrt(2^n n!) (k < 1.0865) for each Cartesian factor:
//      G_m(r) = (2/sqrt(pi)) k^3 sqrt(2^m m!) int_sqrt(p)^inf s^m exp(-s^2 r^2 / 2) ds.
//    It bounds the erfc term itself on the cluster (R >= d - rho) and the
//    terms of its Taylor polynomial at b (R = d).
namespace {
constexpr double cramer = 1.0865;

// G_m(r) for m = 0, ..., values.size() - 1, with the integral written as
// (1/2) (2 / r^2)^((m+1)/2) Gamma((m+1)/2, p r^2 / 2) and the upper
// incomplete gamma function from Gamma(1/2, x) = sqrt(pi) erfc(sqrt(x)),
// Gamma(1, x) = exp(-x) and Gamma(s + 1, x) = s Gamma(s, x) + x^s exp(-x).
void erfc_derivative_bounds(double exponent, double r, std::vector<double>& values) {
    constexpr double pi = std::numbers::pi;
    if (!(r > 0.0)) {
        std::fill(values.begin(), values.end(), std::numeric_limits<double>::infinity());
        return;
    }
    const double x = 0.5 * exponent * r * r;
    const double decay = std::exp(-x);
    const double root_x = std::sqrt(x);
    const double scale = 2.0 / (r * r), root_scale = std::sqrt(scale);
    // Gamma(s, x) for s = 1/2, 3/2, ... (even m) and s = 1, 2, ... (odd m),
    // with x^s and scale^s kept as running products.
    double gamma_half = std::sqrt(pi) * std::erfc(root_x), x_half = root_x;
    double gamma_one = decay, x_one = x;
    double s_half = 0.5, s_one = 1.0;
    double scale_power = root_scale; // scale^((m+1)/2)
    double norm2 = 1.0;              // 2^m m!
    const double prefactor = 2.0 / std::sqrt(pi) * cramer * cramer * cramer * 0.5;
    for (std::size_t m = 0; m < values.size(); ++m) {
        if (m > 0) {
            norm2 *= 2.0 * static_cast<double>(m);
            scale_power *= root_scale;
        }
        double gamma;
        if (m % 2 == 0) {
            if (m > 0) {
                gamma_half = s_half * gamma_half + x_half * decay;
                x_half *= x;
                s_half += 1.0;
            }
            gamma = gamma_half;
        } else {
            if (m > 1) {
                gamma_one = s_one * gamma_one + x_one * decay;
                x_one *= x;
                s_one += 1.0;
            }
            gamma = gamma_one;
        }
        values[m] = prefactor * std::sqrt(norm2) * scale_power * gamma;
    }
}

// T_L(d) for one L; infinite unless the series converges (rho < d).
double taylor_remainder_bound(int L, int K, double rho, double d) {
    // Closer than 1.25 rho the series converges too slowly to be useful;
    // infinity keeps such pairs exact.
    if (!(d > 1.25 * rho)) return std::numeric_limits<double>::infinity();
    const double q = rho / d;
    double term = 1.0 / std::pow(d, L + 1); // n = 0 term without (n+L)!/n!
    for (int j = 1; j <= L; ++j)
        term *= j;
    for (int n = 1; n <= K + 1; ++n)
        term *= (n + L) / static_cast<double>(n) * q;
    double sum = 0.0;
    for (int n = K + 1; n < K + 1000; ++n) {
        sum += term;
        const double ratio = (n + L + 1) / static_cast<double>(n + 1) * q;
        // The ratios decrease with n, so once below 1 the tail is geometric.
        if (ratio < 0.9) return sum + term * ratio / (1.0 - ratio);
        term *= ratio;
        if (!std::isfinite(term)) break;
    }
    return std::numeric_limits<double>::infinity();
}
} // namespace

ElectrostaticPotential::FarFieldPlan
ElectrostaticPotential::far_field_plan(double radius, double tolerance) const {
    if (!std::isfinite(radius) || radius < 0.0)
        throw std::invalid_argument("cluster radius must be finite and nonnegative");
    if (!std::isfinite(tolerance) || tolerance < 0.0)
        throw std::invalid_argument(
            "electrostatic potential tolerance must be finite and nonnegative");
    FarFieldPlan plan;
    plan.radius = radius;
    plan.far_distance2.assign(pairs_.size(), std::numeric_limits<double>::infinity());
    if (tolerance == 0.0 || pairs_.empty()) return plan;
    // The total error at a point is at most the sum of the pair bounds.
    const double allowed = tolerance / static_cast<double>(pairs_.size());
    constexpr int K = local_order;
    std::vector<double> weights, near_bounds, far_bounds;
    for (std::size_t i = 0; i < pairs_.size(); ++i) {
        const auto& pair = pairs_[i];
        const int L = pair.order;
        weights.assign(static_cast<std::size_t>(L + 1), 0.0);
        const double* c = coefficients_.data() + pair.first_coefficient;
        const double to_multipole = 0.5 * std::sqrt(std::numbers::pi / pair.exponent);
        std::size_t k = 0;
        for (int order = 0; order <= L; ++order)
            for (int j = 0; j < (order + 1) * (order + 2) / 2; ++j, ++k)
                weights[static_cast<std::size_t>(order)] += std::abs(c[k]) * to_multipole;
        near_bounds.assign(static_cast<std::size_t>(L + 1), 0.0);
        far_bounds.assign(static_cast<std::size_t>(L + K + 1), 0.0);
        const auto bound = [&](double d) {
            erfc_derivative_bounds(pair.exponent, d - radius, near_bounds);
            erfc_derivative_bounds(pair.exponent, d, far_bounds);
            double total = 0.0;
            for (int order = 0; order <= L; ++order) {
                const double w = weights[static_cast<std::size_t>(order)];
                if (w == 0.0) continue;
                double term = taylor_remainder_bound(order, K, radius, d) +
                              near_bounds[static_cast<std::size_t>(order)];
                double scale = 1.0;
                for (int n = 0; n <= K; ++n) {
                    if (n > 0) scale *= std::sqrt(3.0) * radius / n;
                    term += scale * far_bounds[static_cast<std::size_t>(order + n)];
                }
                total += w * term;
            }
            return total;
        };
        // Smallest distance with bound <= allowed, by bisection; the bound
        // decreases with d. Pairs that never qualify stay exact.
        double high = std::max(2.0 * radius, radius + 1.0 / std::sqrt(pair.exponent));
        while (high < 1e4 && !(bound(high) <= allowed))
            high *= 2.0;
        if (!(bound(high) <= allowed)) continue;
        double low = radius;
        for (int iteration = 0; iteration < 30 && high - low > 1e-2 * high; ++iteration) {
            const double middle = 0.5 * (low + high);
            (bound(middle) <= allowed ? high : low) = middle;
        }
        plan.far_distance2[i] = high * high;
    }
    return plan;
}

void ElectrostaticPotential::evaluate_cluster(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    const Eigen::Vector3d& center, const FarFieldPlan& plan, std::span<double> output) const {
    if (output.size() != static_cast<std::size_t>(points.cols()))
        throw std::invalid_argument(
            "electrostatic potential output size does not match the points");
    if (plan.far_distance2.size() != pairs_.size())
        throw std::invalid_argument("far-field plan does not match the potential");
    if (!center.allFinite()) throw std::invalid_argument("cluster centre must be finite");
    for (Eigen::Index j = 0; j < points.cols(); ++j) {
        if (!points.col(j).allFinite())
            throw std::invalid_argument("electrostatic potential point must be finite");
        // Allow for rounding in the coordinates of the points and centre.
        if ((points.col(j) - center).norm() > plan.radius + 1e-12 * (plan.radius + center.norm()))
            throw std::invalid_argument("cluster point outside the planned radius");
    }
    Scratch scratch;
    prepare(scratch);
    // Far pairs: one local Taylor expansion about the centre,
    //   sum_k c_k R_alpha(P - b - h) = sum_beta h^beta (-1)^|beta| / beta!
    //                                   sum_k c_k R_(alpha+beta)(P - b).
    std::fill(scratch.local.begin(), scratch.local.end(), 0.0);
    std::vector<std::size_t> near;
    for (std::size_t i = 0; i < pairs_.size(); ++i) {
        const auto& pair = pairs_[i];
        if ((pair.center - center).squaredNorm() < plan.far_distance2[i]) {
            near.push_back(i);
            continue;
        }
        hermite_integrals(pair.exponent, pair.center - center, pair.order + local_order, scratch);
        const auto& indices = hermite_indices_[static_cast<std::size_t>(pair.order)];
        const double* c = coefficients_.data() + pair.first_coefficient;
        for (std::size_t k = 0; k < indices.size(); ++k) {
            const double* r = scratch.current.data() + indices[k];
            for (std::size_t b = 0; b < local_offsets_.size(); ++b)
                scratch.local[b] += c[k] * r[local_offsets_[b]];
        }
    }
    for (std::size_t b = 0; b < local_offsets_.size(); ++b)
        scratch.local[b] *= local_factors_[b];

    for (Eigen::Index j = 0; j < points.cols(); ++j) {
        const Eigen::Vector3d point = points.col(j);
        const Eigen::Vector3d h = point - center;
        // Powers of the components of h, then the polynomial.
        std::array<std::array<double, local_order + 1>, 3> powers{};
        for (int a = 0; a < 3; ++a) {
            powers[static_cast<std::size_t>(a)][0] = 1.0;
            for (int e = 1; e <= local_order; ++e)
                powers[static_cast<std::size_t>(a)][static_cast<std::size_t>(e)] =
                    powers[static_cast<std::size_t>(a)][static_cast<std::size_t>(e - 1)] * h[a];
        }
        double polynomial = 0.0;
        for (std::size_t b = 0; b < local_powers_.size(); ++b) {
            const auto& power = local_powers_[b];
            polynomial += scratch.local[b] * powers[0][static_cast<std::size_t>(power[0])] *
                          powers[1][static_cast<std::size_t>(power[1])] *
                          powers[2][static_cast<std::size_t>(power[2])];
        }
        double sum = polynomial;
        for (const std::size_t i : near)
            sum += pair_potential(pairs_[i], point, scratch);
        output[static_cast<std::size_t>(j)] = nuclear(point) - sum;
    }
}
} // namespace chem::kumo::detail
