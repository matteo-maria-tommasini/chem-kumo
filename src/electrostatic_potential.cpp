#include "electrostatic_potential.hpp"
#include "basis_math.hpp"
#include "boys_function.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include <algorithm>
#include <cmath>
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
    const int side = max_order_ + 1;
    hermite_indices_.resize(static_cast<std::size_t>(max_order_ + 1));
    std::vector<std::vector<HermiteIndex>> triples(static_cast<std::size_t>(max_order_ + 1));
    for (int L = 0; L <= max_order_; ++L) {
        triples[static_cast<std::size_t>(L)] = hermite_triples(L);
        for (const auto& h : triples[static_cast<std::size_t>(L)])
            hermite_indices_[static_cast<std::size_t>(L)].push_back((h.t * side + h.u) * side +
                                                                    h.v);
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

// Hermite Coulomb integrals (Helgaker, Jorgensen and Olsen, Molecular
// Electronic-Structure Theory, section 9.9):
//   R^n_000 = (-2p)^n F_n(p |PC|^2),
//   R^n_(t+1,u,v) = t R^(n+1)_(t-1,u,v) + X_PC R^(n+1)_(t,u,v),
// and likewise for u and v, from n = L down to n = 0. The integral of a
// Hermite Gaussian Lambda_tuv over 1/|r - C| is (2 pi / p) R^0_tuv; the factor
// 2 pi / p is already in the stored coefficients.
double ElectrostaticPotential::electronic(const Eigen::Vector3d& point, Scratch& scratch) const {
    if (!point.allFinite())
        throw std::invalid_argument("electrostatic potential point must be finite");
    const int side = max_order_ + 1;
    const auto cube = static_cast<std::size_t>(side * side * side);
    if (scratch.current.size() != cube) {
        scratch.boys.assign(static_cast<std::size_t>(side), 0.0);
        scratch.previous.assign(cube, 0.0);
        scratch.current.assign(cube, 0.0);
        scratch.triples = hermite_triples(max_order_);
    }
    const auto offset = [side](const HermiteIndex& h) {
        return static_cast<std::size_t>((h.t * side + h.u) * side + h.v);
    };
    double sum = 0.0;
    for (const auto& pair : pairs_) {
        const Eigen::Vector3d pc = pair.center - point;
        const int L = pair.order;
        const std::span<double> boys(scratch.boys.data(), static_cast<std::size_t>(L + 1));
        boys_function(pair.exponent * pc.squaredNorm(), boys);
        // Triples with t + u + v <= s are the first C(s + 3, 3) entries.
        const auto count = [](int s) {
            return static_cast<std::size_t>((s + 1) * (s + 2) * (s + 3) / 6);
        };
        double power = 1.0;
        for (int n = 0; n < L; ++n)
            power *= -2.0 * pair.exponent;
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
                    const HermiteIndex lower{h.t - 1, h.u, h.v};
                    value = pc.x() * previous[offset(lower)];
                    if (h.t > 1) value += (h.t - 1) * previous[offset({h.t - 2, h.u, h.v})];
                } else if (h.u > 0) {
                    const HermiteIndex lower{h.t, h.u - 1, h.v};
                    value = pc.y() * previous[offset(lower)];
                    if (h.u > 1) value += (h.u - 1) * previous[offset({h.t, h.u - 2, h.v})];
                } else {
                    const HermiteIndex lower{h.t, h.u, h.v - 1};
                    value = pc.z() * previous[offset(lower)];
                    if (h.v > 1) value += (h.v - 1) * previous[offset({h.t, h.u, h.v - 2})];
                }
                current[offset(h)] = value;
            }
            if (n > 0) power /= -2.0 * pair.exponent;
        }
        const auto& indices = hermite_indices_[static_cast<std::size_t>(L)];
        const double* c = coefficients_.data() + pair.first_coefficient;
        for (std::size_t k = 0; k < indices.size(); ++k)
            sum += c[k] * scratch.current[static_cast<std::size_t>(indices[k])];
    }
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
} // namespace chem::kumo::detail
