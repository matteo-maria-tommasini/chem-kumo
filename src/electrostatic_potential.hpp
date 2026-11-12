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
///   V(r) = sum_A Z_A / |r - R_A| - sum_(mu,nu) P_(mu nu) integral chi_mu chi_nu / |r - r'| dr',
/// with the nuclear charges Z_A of the atoms and the AO density matrix P in
/// Gaussian FCHK AO order.
///
/// The electronic term uses the McMurchie-Davidson scheme. Each product of
/// two Cartesian primitives is a sum of Hermite Gaussians centred at the
/// Gaussian product centre; P is contracted with their expansion
/// coefficients once, at construction, so that each point needs, per
/// primitive pair, one Boys function evaluation, the Hermite Coulomb
/// integrals R_tuv and a dot product. Pure shells are handled by
/// transforming P to Cartesian components. The pointwise functions
/// (operator(), electronic(), evaluate()) include every primitive pair, so
/// they give the exact integral up to rounding; evaluate_cluster() replaces
/// pairs far from a cluster of points by a local Taylor expansion, within a
/// tolerance set by its far-field plan.
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
    /// finite, or a shell has an angular momentum above 16.
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

    /// Precomputed data for clusters of points within `radius` of a centre.
    struct FarFieldPlan {
        static constexpr int level_count = 33;
        static constexpr double cell_size = 2.0; ///< bohr, edge of the source cells
        double radius = 0.0;
        double tolerance = 0.0;
        /// Budget of the individually expanded pairs and its error levels,
        /// pair_tolerance * 10^(-j/2).
        double pair_tolerance = 0.0;
        std::array<double, level_count> levels{};
        /// Per pair and level: squared distance from a cluster centre beyond
        /// which the expansion error of the pair is at most that level
        /// (infinite if never).
        std::vector<float> level_distance2;

        /// Pairs grouped by cell and exponent, with their multipole moments.
        struct SourceCluster {
            Eigen::Vector3d center = Eigen::Vector3d::Zero();
            double radius = 0.0;       ///< largest distance of a pair centre
            double min_distance = 0.0; ///< translate only from at least this far
            std::vector<std::size_t> pairs;
            std::vector<double> moments; ///< (multipole_order + 1) x moment_count
        };

        std::vector<SourceCluster> clusters;
        std::size_t moment_count = 0;
    };

    /// Plan for clusters of points within `radius` (bohr) of their centre,
    /// with an error of at most `tolerance` (hartree/e) at every point; a
    /// tolerance of 0 keeps every pair exact. See electrostatic_potential.cpp
    /// for the error bounds.
    [[nodiscard]] FarFieldPlan far_field_plan(double radius, double tolerance) const;

    /// Total potential at the columns of `points`, which must lie within
    /// plan.radius of `center`. Source clusters far from the centre contribute
    /// through one multipole-to-local translation each; of the other pairs,
    /// those whose error bounds fit in the budget through individual Taylor
    /// expansions about the centre; the rest exactly. Each value depends only
    /// on the point, the centre and the plan.
    void evaluate_cluster(const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                          const Eigen::Vector3d& center, const FarFieldPlan& plan,
                          std::span<double> output) const;

private:
    /// Hermite expansion of one primitive pair, already contracted with P:
    /// the electronic potential is -sum_tuv coefficients[k] R_tuv(p, P - r),
    /// with (t, u, v) the k-th entry of hermite_indices_[order].
    struct PrimitivePair {
        double exponent = 0.0;             ///< p = alpha + beta
        Eigen::Vector3d center;            ///< Gaussian product centre
        int order = 0;                     ///< la + lb
        std::size_t first_coefficient = 0; ///< into coefficients_
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

    std::vector<PrimitivePair> pairs_;
    std::vector<double> coefficients_;
    /// For each order L, the offsets of (t, u, v), t + u + v <= L, in a
    /// side_^3 table.
    std::vector<std::vector<int>> hermite_indices_;
    int max_order_ = 0; ///< largest la + lb
    int side_ = 1;      ///< max_order_ + multipole_order + 1
    std::vector<int> local_offsets_;
    std::vector<int> moment_offsets_;
    std::vector<double> local_factors_;
    std::vector<std::array<int, 3>> local_powers_;
    std::vector<Atom> atoms_;
};
} // namespace chem::kumo::detail
