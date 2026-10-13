#include "chem/kumo.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Points = Eigen::Matrix<double, 3, Eigen::Dynamic>;

double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

// Three-dimensional lattice with the same point count and coordinates for all
// HBC basis sets. The bounding box depends only on the nuclear geometry.
struct SampleGrid {
    chem::kumo::Grid grid;
    Points points;
};

SampleGrid make_grid(const std::vector<chem::kumo::Atom>& atoms,
                     std::size_t requested, double margin) {
    if (atoms.empty()) throw std::invalid_argument("empty molecule");
    const double side = std::ceil(std::cbrt(static_cast<double>(requested)));
    if (side > 1000.0) throw std::invalid_argument("too many points");
    const int n = static_cast<int>(side);
    chem::kumo::Grid grid;
    grid.counts = Eigen::Vector3i(n, n, n);
    Eigen::Vector3d lo = atoms.front().position;
    Eigen::Vector3d hi = lo;
    for (const auto& atom : atoms) {
        lo = lo.cwiseMin(atom.position);
        hi = hi.cwiseMax(atom.position);
    }
    lo.array() -= margin;
    hi.array() += margin;
    grid.origin = lo;
    grid.steps.setZero();
    for (int k = 0; k < 3; ++k)
        grid.steps(k,k) = (hi[k] - lo[k]) / (n - 1);
    const auto count = static_cast<Eigen::Index>(n) * n * n;
    Points points(3, count);
    Eigen::Index q = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                points.col(q++) = grid.origin + grid.steps * Eigen::Vector3d(i,j,k);
    return {grid, std::move(points)};
}

// Consume returned values so the compiler cannot discard the benchmark work.
// The checksum is computed outside timed regions.
double sample(const std::vector<chem::kumo::DensityDerivatives>& data) {
    const auto& d = data.at(data.size()/2);
    return d.value + d.gradient.squaredNorm() + d.hessian.squaredNorm()
           + d.laplacian + d.tau;
}

// Report minimum wall time across repeats, with one warm-up. No timing-based
// pass/fail conditions: this is a developer benchmark, not a CI test.
template <typename F>
double measure(int repeats, F&& f) {
    f(); // warm-up, excluded
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < repeats; ++i) {
        const auto start = Clock::now();
        f();
        best = std::min(best, elapsed(start));
    }
    return best;
}
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 5) {
        std::cerr << "Usage: kumo-spatial-benchmark FILE.fchk[.gz] [points=125] [repeats=2] [margin_bohr=2.5]\n";
        return 2;
    }
    try {
        const std::size_t requested = argc > 2 ? std::stoull(argv[2]) : 125;
        const int repeats = argc > 3 ? std::stoi(argv[3]) : 2;
        const double margin = argc > 4 ? std::stod(argv[4]) : 2.5;
        if (requested < 8 || requested > 1000000 || repeats < 1 || repeats > 20 ||
            !std::isfinite(margin) || margin <= 0.0)
            throw std::invalid_argument("invalid benchmark arguments");

        const auto mol = chem::kumo::Molecule::from_fchk(argv[1]);
        const auto data = make_grid(mol.atoms(), requested, margin);
        const auto count = static_cast<std::size_t>(data.points.cols());
        const chem::kumo::FormattedCheckpoint fchk(argv[1]);
        const chem::kumo::BasisEvaluator evaluator{chem::kumo::BasisSet(fchk)};
        double checksum = 0.0;

        // AO timings exclude parsing, construction, and grid setup.
        Eigen::VectorXd ao;
        const double ao_values = measure(repeats, [&] {
            for (Eigen::Index i = 0; i < data.points.cols(); ++i)
                ao = evaluator.evaluate(data.points.col(i));
        });
        checksum += ao.squaredNorm();

        chem::kumo::AODerivatives deriv;
        const double ao_derivatives = measure(repeats, [&] {
            for (Eigen::Index i = 0; i < data.points.cols(); ++i)
                deriv = evaluator.evaluate_derivatives(data.points.col(i));
        });
        checksum += deriv.values.squaredNorm() + deriv.gradient.squaredNorm();

        // Isolate the four dense matrix products required by the derivative
        // contractions, using AO data prepared OUTSIDE the timed section.
        // This is a representative block, not an estimate of total runtime.
        const Eigen::Index b = std::min<Eigen::Index>(64, data.points.cols());
        const Eigen::Index nbas = static_cast<Eigen::Index>(mol.basis_size());
        Eigen::MatrixXd phi(nbas,b), gx(nbas,b), gy(nbas,b), gz(nbas,b);
        for (Eigen::Index j = 0; j < b; ++j) {
            const auto d = evaluator.evaluate_derivatives(data.points.col(j));
            phi.col(j) = d.values;
            gx.col(j) = d.gradient.col(0);
            gy.col(j) = d.gradient.col(1);
            gz.col(j) = d.gradient.col(2);
        }
        const auto& packed = fchk.reals("Total SCF Density");
        if (packed.size() != static_cast<std::size_t>(nbas*(nbas+1)/2))
            throw std::runtime_error("invalid packed density size");
        Eigen::MatrixXd p(nbas,nbas);
        std::size_t offset = 0;
        for (Eigen::Index i = 0; i < nbas; ++i)
            for (Eigen::Index j = 0; j <= i; ++j)
                p(i,j) = p(j,i) = packed[offset++];
        Eigen::MatrixXd yp(nbas,b), yx(nbas,b), yy(nbas,b), yz(nbas,b);
        const double contraction = measure(repeats, [&] {
            yp.noalias() = p*phi;
            yx.noalias() = p*gx;
            yy.noalias() = p*gy;
            yz.noalias() = p*gz;
        });
        checksum += yp(0,0) + yx(0,0) + yy(0,0) + yz(0,0);

        // A/B comparison: same AO panels and density matrix, but a single
        // Eigen matrix multiplication instead of four separate products.
        Eigen::MatrixXd packed_ao(nbas,4*b), fused(nbas,4*b);
        packed_ao.middleCols(0,b) = phi;
        packed_ao.middleCols(b,b) = gx;
        packed_ao.middleCols(2*b,b) = gy;
        packed_ao.middleCols(3*b,b) = gz;
        const double fused_contraction = measure(repeats, [&] {
            fused.noalias() = p * packed_ao;
        });
        const double fused_max_abs = std::max({
            (fused.middleCols(0,b)-yp).cwiseAbs().maxCoeff(),
            (fused.middleCols(b,b)-yx).cwiseAbs().maxCoeff(),
            (fused.middleCols(2*b,b)-yy).cwiseAbs().maxCoeff(),
            (fused.middleCols(3*b,b)-yz).cwiseAbs().maxCoeff()});
        if (fused_max_abs > 1.e-9)
            throw std::runtime_error("fused and separate contractions disagree");
        checksum += fused(0,0);

        // Density-only reference uses the existing blocked Eigen implementation.
        chem::kumo::Field density = mol.density(data.grid, chem::kumo::Spin::Total);
        const double density_only = measure(repeats, [&] {
            density = mol.density(data.grid, chem::kumo::Spin::Total);
        });
        checksum += density.values()[count/2];

        chem::kumo::DensityDerivatives pointwise;
        const double pointwise_time = measure(repeats, [&] {
            for (Eigen::Index i = 0; i < data.points.cols(); ++i)
                pointwise = mol.density_derivatives(data.points.col(i), chem::kumo::Spin::Total);
        });
        checksum += pointwise.value + pointwise.tau;

        std::cout << std::setprecision(12)
                  << "file=" << argv[1] << "\n"
                  << "basis=" << mol.basis_size() << " atoms=" << mol.atoms().size()
                  << " points=" << count << " repeats=" << repeats
                  << " margin_bohr=" << margin << "\n"
                  << "ao_values_s=" << ao_values << "\n"
                  << "ao_derivatives_s=" << ao_derivatives << "\n"
                  << "ao_derivatives_over_values=" << ao_derivatives/ao_values << "\n"
                  << "contraction_block_points=" << b << "\n"
                  << "four_dense_contractions_s=" << contraction << "\n"
                  << "fused_dense_contraction_s=" << fused_contraction << "\n"
                  << "fused_contraction_speedup=" << contraction/fused_contraction << "\n"
                  << "fused_contraction_max_abs=" << fused_max_abs << "\n"
                  << "density_only_s=" << density_only << "\n"
                  << "derivatives_pointwise_s=" << pointwise_time << "\n"
                  << "derivatives_over_density=" << pointwise_time/density_only << "\n";

        for (std::size_t block : {16u, 32u, 64u, 128u, 256u}) {
            std::vector<chem::kumo::DensityDerivatives> result;
            const double t = measure(repeats, [&] {
                result = mol.density_derivatives_blocked(
                    data.points, chem::kumo::Spin::Total, block);
            });
            checksum += sample(result);
            std::cout << "derivatives_block_" << block << "_s=" << t
                      << " speedup=" << pointwise_time/t << "\n";
        }
        std::cout << "checksum=" << checksum << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
