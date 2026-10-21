#include "chem/kumo/local_density_overlap_multichannel_hierarchical.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto mol = Molecule::from_fchk(argv[arg]);
            const Eigen::MatrixXd total = mol.density_matrix(Spin::Total);
            const Eigen::MatrixXd spin = mol.density_matrix(Spin::Difference);
            const Eigen::MatrixXd alpha = 0.5 * (total + spin);
            const Eigen::MatrixXd beta = 0.5 * (total - spin);
            const Eigen::Index n = total.rows();
            const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(n, n);
            Eigen::MatrixXd artificial_a(n,n), artificial_b(n,n);
            for (Eigen::Index i=0;i<n;++i)
                for (Eigen::Index j=0;j<n;++j) {
                    artificial_a(i,j) = 0.15*std::cos(0.17*double(i+j+2)) + (i==j ? 0.4 : 0.0);
                    artificial_b(i,j) = 0.11*std::sin(0.13*double((i+1)*(j+1))) + (i==j ? 0.2 : 0.0);
                }
            const std::vector<std::vector<Eigen::MatrixXd>> cases{
                {alpha,beta,zero}, {artificial_a,artificial_b,zero}, {zero,zero}
            };
            const std::vector<std::vector<double>> budgets{
                {1e-8,1e-10,0.0}, {1e-8,1e-10,0.0}, {0.0,0.0}
            };
            for (std::size_t c=0;c<cases.size();++c) {
                const auto exact = local_density_overlap_multichannel(basis,cases[c]);
                const auto result = local_density_overlap_multichannel_hierarchical(
                    basis,cases[c],budgets[c],1024*1024);
                const auto& s=result.statistics;
                if (s.evaluated_quartets+s.skipped_quartets!=s.total_quartets ||
                    s.candidate_visits!=s.evaluated_quartets)
                    throw std::runtime_error("hierarchical quartet accounting mismatch");
                double worst=0;
                for (std::size_t k=0;k<cases[c].size();++k) {
                    const double scale=std::max(1.0,exact.matrices[k].cwiseAbs().maxCoeff());
                    const double error=(result.matrices[k]-exact.matrices[k]).cwiseAbs().maxCoeff();
                    const double bound=s.max_error_bounds[k];
                    if (!std::isfinite(bound) || bound>budgets[c][k]+1e-18 ||
                        error>budgets[c][k]+2e-11*scale ||
                        (result.matrices[k]-result.matrices[k].transpose()).cwiseAbs().maxCoeff()>1e-11*scale)
                        throw std::runtime_error("hierarchical accuracy or budget failure");
                    worst=std::max(worst,error);
                }
                // Compare the physical-case quartet count with the independent
                // bidirectional cutoff diagnostic (same alpha/beta budgets).
                if (c==0 && n==28 && s.skipped_quartets!=7)
                    throw std::runtime_error("O2 cutoff count mismatch");
                if (c==0 && n==64 && s.skipped_quartets!=11146)
                    throw std::runtime_error("acrolein cutoff count mismatch");
                std::cout << "system=" << argv[arg] << " case=" << c
                          << " skipped=" << s.skipped_quartets
                          << " evaluated=" << s.evaluated_quartets
                          << " candidate_visits=" << s.candidate_visits
                          << " binary_nodes=" << s.binary_nodes
                          << " max_abs=" << worst
                          << " total_seconds=" << s.total_seconds << '\n';
            }
        }
    } catch(const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n'; return 1;
    }
}
