#include "chem/kumo/local_density_overlap_multichannel_hierarchical.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

template<class Fn>
static void rejects(Fn&& fn, const char* message) {
    bool rejected = false;
    try { fn(); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, message);
}

static void validate(const LocalHierarchicalDensityOverlap& result,
                     const LocalMultichannelDensityOverlap& reference,
                     const std::vector<double>& eps,
                     const std::string& label) {
    const auto& s = result.statistics;
    require(s.evaluated_quartets + s.skipped_quartets == s.total_quartets,
            "hierarchical quartet accounting");
    require(s.candidate_visits == s.evaluated_quartets, "hierarchical candidate accounting");
    require(result.matrices.size() == eps.size(), "hierarchical output count");
    for (std::size_t k = 0; k < eps.size(); ++k) {
        const auto& w = result.matrices[k];
        const auto& r = reference.matrices[k];
        const double scale = std::max(1.0, r.cwiseAbs().maxCoeff());
        const double err = (w-r).cwiseAbs().maxCoeff();
        const double sym = (w-w.transpose()).cwiseAbs().maxCoeff();
        const double bound = s.max_error_bounds.at(k);
        require(w.allFinite() && std::isfinite(bound), "nonfinite hierarchical result");
        require(bound >= 0 && bound <= eps[k] + 1e-18, "hierarchical budget bound");
        require(err <= eps[k] + 2e-11*scale, "hierarchical numerical error");
        require(sym <= 1e-11*scale, "hierarchical symmetry");
        std::cout << " case=" << label << " channel=" << k
                  << " epsilon=" << eps[k] << " max_abs=" << err
                  << " max_bound=" << bound;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        const bool benchmark = std::string(argv[1]) == "--benchmark";
        const int first = benchmark ? 2 : 1;
        if (first >= argc) return 2;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = first; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto mol = Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd total = mol.density_matrix(Spin::Total);
            const Eigen::MatrixXd spin = mol.density_matrix(Spin::Difference);
            const Eigen::MatrixXd alpha = 0.5 * (total + spin);
            const Eigen::MatrixXd beta = 0.5 * (total - spin);
            const Eigen::Index n = total.rows();
            const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(n,n);
            Eigen::MatrixXd artificial(n,n);
            for (Eigen::Index i=0;i<n;++i)
                for (Eigen::Index j=0;j<n;++j)
                    artificial(i,j) = 0.13*std::cos(0.19*double(i+j+2)) + (i==j ? 0.35 : 0.0);
            const std::vector<Eigen::MatrixXd> physical{alpha,beta,zero};
            const std::vector<double> eps{1e-8,1e-10,0.0};
            constexpr std::size_t cap = 1024*1024;
            const auto reference = local_density_overlap_multichannel(basis,physical);
            const auto hierarchical = local_density_overlap_multichannel_hierarchical(basis,physical,eps,cap);
            validate(hierarchical,reference,eps,"physical");
            std::cout << " system=" << argv[arg] << " skipped="
                      << hierarchical.statistics.skipped_quartets << " OK\n";

            // Exact zero budgets and duplicate channels exercise independent
            // per-channel validation and the zero-screening boundary.
            const std::vector<Eigen::MatrixXd> duplicated{artificial,artificial,zero};
            const std::vector<double> strict{0.0,1e-8,0.0};
            const auto exact_duplicate = local_density_overlap_multichannel(basis,duplicated);
            const auto strict_result = local_density_overlap_multichannel_hierarchical(
                basis,duplicated,strict,cap);
            validate(strict_result,exact_duplicate,strict,"zero_budget_duplicate");
            std::cout << " system=" << argv[arg] << " strict=OK\n";

            const std::vector<Eigen::MatrixXd> all_zero{zero,zero};
            const auto zeros = local_density_overlap_multichannel_hierarchical(
                basis,all_zero,{0.0,0.0},cap);
            require(zeros.statistics.evaluated_quartets == 0 &&
                    zeros.statistics.candidate_visits == 0 &&
                    zeros.statistics.active_channels == 0,
                    "inactive channel fast path");
            require(zeros.matrices[0].isZero(0.0) && zeros.matrices[1].isZero(0.0),
                    "inactive matrices not zero");
            rejects([&]{ (void)local_density_overlap_multichannel_hierarchical(
                basis,physical,{1e-8},cap); }, "budget count accepted");
            rejects([&]{ (void)local_density_overlap_multichannel_hierarchical(
                basis,physical,{-1.0,1e-8,0.0},cap); }, "negative budget accepted");
            rejects([&]{ (void)local_density_overlap_multichannel_hierarchical(
                basis,physical,{std::numeric_limits<double>::quiet_NaN(),1e-8,0.0},cap); },
                "NaN budget accepted");
            Eigen::MatrixXd wrong = Eigen::MatrixXd::Zero(n+1,n+1);
            rejects([&]{ (void)local_density_overlap_multichannel_hierarchical(
                basis,{wrong},{0.0},cap); }, "wrong dimensions accepted");
            Eigen::MatrixXd asymmetric = artificial;
            if (n > 1) asymmetric(0,1) += 0.1;
            rejects([&]{ (void)local_density_overlap_multichannel_hierarchical(
                basis,{asymmetric},{1e-8},cap); }, "asymmetry accepted");
            Eigen::MatrixXd nonfinite = artificial;
            nonfinite(0,0) = std::numeric_limits<double>::quiet_NaN();
            rejects([&]{ (void)local_density_overlap_multichannel_hierarchical(
                basis,{nonfinite},{1e-8},cap); }, "nonfinite density accepted");
            std::cout << " system=" << argv[arg] << " edge_cases=OK\n";

            if (benchmark) {
                // Each backend runs in the same process on the same inputs.
                // Timing is diagnostic only; no performance assertion.
                for (int rep=0;rep<3;++rep) {
                    const auto start = std::chrono::steady_clock::now();
                    const auto canonical = local_density_overlap_multichannel_screened(
                        basis,physical,eps,cap,MultichannelScreeningOrder::Canonical);
                    const auto maximum = local_density_overlap_multichannel_screened(
                        basis,physical,eps,cap,MultichannelScreeningOrder::MaximumImportance);
                    const auto hierarchical_run = local_density_overlap_multichannel_hierarchical(
                        basis,physical,eps,cap);
                    const double wall = std::chrono::duration<double>(
                        std::chrono::steady_clock::now()-start).count();
                    for (std::size_t k=0;k<physical.size();++k) {
                        const double scale=std::max(1.0,reference.matrices[k].cwiseAbs().maxCoeff());
                        for (const auto* w : {&canonical.matrices[k], &maximum.matrices[k],
                                              &hierarchical_run.matrices[k]})
                            require((*w-reference.matrices[k]).cwiseAbs().maxCoeff()
                                    <= eps[k]+2e-11*scale,"benchmark numerical regression");
                    }
                    std::cout << " benchmark=" << argv[arg] << " repetition=" << rep
                              << " canonical_seconds=" << canonical.statistics.total_seconds
                              << " maximum_seconds=" << maximum.statistics.total_seconds
                              << " hierarchical_seconds=" << hierarchical_run.statistics.total_seconds
                              << " canonical_skipped=" << canonical.statistics.skipped_quartets
                              << " maximum_skipped=" << maximum.statistics.skipped_quartets
                              << " hierarchical_skipped=" << hierarchical_run.statistics.skipped_quartets
                              << " combined_wall_seconds=" << wall << '\n';
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
