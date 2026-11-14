#pragma once
#include "chem/kumo/basis.hpp"
#include "chem/kumo/density.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace chem::kumo::detail {
/// Molecular electrostatic potential, in hartree/e (atomic units),
///   V(r) = sum_A Z_A / |r - R_A|
///          - sum_(mu,nu) P_(mu nu) integral chi_mu(r') chi_nu(r') / |r - r'| dr',
/// with the nuclear charges Z_A and positions R_A of the atoms, the atomic
/// orbitals (AOs) chi_mu and the AO density matrix P in Gaussian FCHK AO
/// order. The mathematics is in 'Molecular electrostatic potential' in
/// docs/mathematics/chem-kumo-mathematics.tex.
///
/// The electronic term uses the McMurchie-Davidson scheme. Each product of
/// two Cartesian primitives is a sum of Hermite Gaussians centred at the
/// Gaussian product centre; P is contracted with their expansion
/// coefficients once, at construction, so that each point needs, per
/// primitive pair, one Boys function evaluation, the Hermite Coulomb
/// integrals R_tuv and a dot product. Pure (spherical) shells are handled by
/// transforming P to Cartesian components. The pointwise functions
/// (operator(), electronic(), evaluate()) include every primitive pair, so
/// they give the exact integral up to rounding; evaluate_cluster() replaces
/// pairs far from a cluster of points by Taylor expansions about the cluster
/// centre, within the tolerance of its FarFieldPlan.
///
/// At a point closer than nucleus_exclusion_radius to a nucleus, that
/// nucleus's own term, singular there, is left out: the value is the
/// potential that this nucleus feels from the electrons and the other nuclei
/// (the electrostatic potential at the nucleus). Regular grids often contain
/// nuclear positions exactly, for example kumo's default grids.
class ElectrostaticPotential {
public:
    static constexpr double nucleus_exclusion_radius = 1e-8; ///< bohr

    /// Throws std::invalid_argument if P does not match the basis, is not
    /// finite, or a shell has an angular momentum above 13.
    ElectrostaticPotential(const BasisSet& basis, const Eigen::MatrixXd& density_matrix,
                           std::vector<Atom> atoms);

    /// Electronic term, -integral rho(r') / |r - r'| dr'.
    [[nodiscard]] double electronic(const Eigen::Vector3d& point) const;
    /// Nuclear term, sum_A Z_A / |r - R_A|, without a nucleus at the point.
    [[nodiscard]] double nuclear(const Eigen::Vector3d& point) const;

    /// Total potential, nuclear + electronic.
    [[nodiscard]] double operator()(const Eigen::Vector3d& point) const {
        return nuclear(point) + electronic(point);
    }

    /// Total potential at the columns of `points`, written to `output`.
    void evaluate(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                  std::span<double> output) const;

    [[nodiscard]] std::size_t primitive_pair_count() const noexcept { return pairs_.size(); }

    /// Order K of the local Taylor expansions of individual primitive pairs.
    static constexpr int local_order = 6;
    /// Total order N of the multipole-to-local expansions of source clusters.
    static constexpr int multipole_order = 10;

    /// Precomputed data for evaluate_cluster(), valid for every cluster of
    /// target points within `radius` of its centre. See 'Far-field
    /// approximations' in docs/mathematics/chem-kumo-mathematics.tex.
    struct FarFieldPlan {
        static constexpr int level_count = 33;
        static constexpr double cell_size = 2.0; ///< bohr, edge of the source cells
        double radius = 0.0;                     ///< bohr, radius of the target clusters
        double tolerance = 0.0;                  ///< hartree/e, largest error at any point
        /// Share of the tolerance for the individually expanded pairs.
        double pair_tolerance = 0.0;
        /// Error levels levels[j] = pair_tolerance * 10^(-j/2).
        std::array<double, level_count> levels{};
        /// For pair i and level j, entry i * level_count + j: the squared
        /// distance from a cluster centre beyond which the expansion error of
        /// the pair is at most levels[j] (infinite if never).
        std::vector<float> level_distance2;

        /// Primitive pairs whose centres lie in the same cube of edge
        /// cell_size and whose exponents lie in the same interval
        /// [4^k, 4^(k+1)), with their multipole moments about `center`.
        struct SourceCluster {
            Eigen::Vector3d center = Eigen::Vector3d::Zero();
            double radius = 0.0; ///< largest distance of a pair centre from `center`
            /// The multipole expansion is used only for target clusters whose
            /// centre is at least this far from `center`; nearer ones treat
            /// the pairs one by one.
            double min_distance = 0.0;
            std::vector<std::size_t> pairs; ///< indices into the primitive pairs
            /// moment_count moments for each m = 0, ..., multipole_order; the
            /// block of m holds the moments truncated at order m.
            std::vector<double> moments;
        };

        std::vector<SourceCluster> clusters;
        std::size_t moment_count = 0;
    };

    /// Plan for clusters of points within `radius` (bohr) of their centre,
    /// with an error of at most `tolerance` (hartree/e) at every point; a
    /// tolerance of 0 keeps every pair exact. The error bounds are derived in
    /// 'Far-field approximations' in docs/mathematics/chem-kumo-mathematics.tex.
    [[nodiscard]] FarFieldPlan far_field_plan(double radius, double tolerance) const;

    /// Total potential at the columns of `points`, which must lie within
    /// plan.radius of `center`. Each source cluster far enough from the centre
    /// contributes through one multipole-to-local translation. Each remaining
    /// pair is replaced by its own Taylor expansion about the centre if its
    /// error bound still fits in the error budget, and is computed exactly
    /// otherwise. Each value depends only on the point, the centre and the
    /// plan, not on the other points.
    void evaluate_cluster(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                          const Eigen::Vector3d& center, const FarFieldPlan& plan,
                          std::span<double> output) const;

private:
    /// Hermite expansion of one pair of primitive Gaussians (exponents alpha
    /// and beta, angular momenta la and lb), already contracted with P. Its
    /// contribution to the electronic potential at r is
    /// -sum_k coefficients_[first_coefficient + k] R_tuv(p, Q - r), with Q the
    /// pair centre and (t, u, v) the k-th entry of hermite_indices_[order].
    struct PrimitivePair {
        double exponent = 0.0;             ///< p = alpha + beta
        Eigen::Vector3d center;            ///< Q, the Gaussian product centre
        int order = 0;                     ///< la + lb
        std::size_t first_coefficient = 0; ///< index of the first coefficient in coefficients_
    };

    struct Scratch;
    void prepare(Scratch& scratch) const;
    void hermite_integrals(double exponent, const Eigen::Vector3d& pc, int L,
                           Scratch& scratch) const;
    void coulomb_derivatives(const Eigen::Vector3d& x, int L, Scratch& scratch) const;
    void hermite_recursion(const Eigen::Vector3d& pc, int L, Scratch& scratch) const;
    [[nodiscard]] double pair_potential(const PrimitivePair& pair, const Eigen::Vector3d& point,
                                        Scratch& scratch) const;
    [[nodiscard]] double electronic(const Eigen::Vector3d& point, Scratch& scratch) const;
    /// weights[L] = W_L, the sum of |e_k| over the Hermite terms of order
    /// L = t + u + v of a pair, with e_k its point-multipole weights.
    void multipole_weights(const PrimitivePair& pair, std::vector<double>& weights) const;
    void plan_pair_levels(FarFieldPlan& plan) const;
    void plan_source_clusters(FarFieldPlan& plan) const;
    [[nodiscard]] std::vector<std::size_t> expand_far_field(const Eigen::Vector3d& center,
                                                            const FarFieldPlan& plan,
                                                            Scratch& scratch) const;
    [[nodiscard]] double local_polynomial(const Eigen::Vector3d& h, const Scratch& scratch) const;
    /// Offset of (t, u, v) in a cubic table of side_ entries per index.
    [[nodiscard]] std::size_t hermite_offset(int t, int u, int v) const;

    std::vector<PrimitivePair> pairs_;
    std::vector<double> coefficients_; ///< Hermite coefficients of all pairs
    /// For each order L, the table offsets (hermite_offset()) of the
    /// (t, u, v) with t + u + v <= L, ordered by t + u + v.
    std::vector<std::vector<int>> hermite_indices_;
    int max_order_ = 0; ///< largest la + lb
    int side_ = 1;      ///< max_order_ + multipole_order + 1
    /// Terms beta = (t, u, v), |beta| <= multipole_order, of the local
    /// expansion: their table offsets, factors (-1)^|beta| / beta! and powers.
    std::vector<int> local_offsets_;
    /// Table offsets of the moment indices nu, |nu| <= max_order_ + multipole_order.
    std::vector<int> moment_offsets_;
    std::vector<double> local_factors_;
    std::vector<std::array<int, 3>> local_powers_;
    std::vector<Atom> atoms_;
};
} // namespace chem::kumo::detail
