#include "chem/kumo/local_density_overlap.hpp"
#include "chem/kumo/local_density_overlap_symmetric.hpp"
#include "chem/kumo/local_density_overlap_symmetric_paired.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace chem::kumo;

namespace {
double max_abs(const Eigen::MatrixXd& a) {
    return a.size() ? a.cwiseAbs().maxCoeff() : 0.0;
}

void compare(const BasisSet& basis, const Eigen::MatrixXd& p,
             const std::string& system, const std::string& name) {
    const Eigen::MatrixXd ref = local_density_overlap(basis, p);
    const Eigen::MatrixXd sym = local_density_overlap_symmetric(basis, p);
    const Eigen::MatrixXd paired = local_density_overlap_symmetric_paired(basis, p);
    const double scale = std::max(1.0, max_abs(ref));
    const double e_sym = max_abs(sym - ref);
    const double e_paired = max_abs(paired - ref);
    const double tolerance = 1e-11 * scale;
    std::cout << std::scientific << std::setprecision(12)
              << "system=" << system << " case=" << name
              << " nao=" << p.rows() << " scale=" << scale
              << " symmetric_error=" << e_sym
              << " paired_error=" << e_paired
              << " tolerance=" << tolerance << '\n';
    if (!ref.allFinite() || !sym.allFinite() || !paired.allFinite() ||
        e_sym > tolerance || e_paired > tolerance)
        throw std::runtime_error("AO pair symmetry regression: " + system + " / " + name);
}
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: test <fchk.gz> [fchk.gz ...]\n";
        return 2;
    }
    try {
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const auto shells = local_shell_components(basis);
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            for (std::size_t a = 0; a < shells.size(); ++a)
                offsets[a + 1] = offsets[a] + static_cast<Eigen::Index>(shells[a].size());
            const Eigen::Index nao = offsets.back();
            if (nao < 2) throw std::runtime_error("need at least two AOs");
            std::size_t shell = shells.size();
            for (std::size_t a = 0; a < shells.size(); ++a)
                if (shells[a].size() >= 2) { shell = a; break; }
            if (shell == shells.size())
                throw std::runtime_error("fixture has no multi-component shell");
            const Eigen::Index x = offsets[shell];
            const Eigen::Index y = x + 1;
            Eigen::MatrixXd p = Eigen::MatrixXd::Zero(nao, nao);
            p(x, y) = p(y, x) = 0.375;
            compare(basis, p, argv[arg], "same_shell_offdiagonal");

            p.setZero();
            p(x, x) = 0.75;
            compare(basis, p, argv[arg], "single_diagonal");

            p.setZero();
            for (Eigen::Index i = 0; i < nao; ++i) {
                p(i, i) = 0.3 + 0.01 * static_cast<double>(i % 7);
                for (Eigen::Index j = i + 1; j < nao; ++j) {
                    const double v = 0.07 * std::sin(0.31 * double(i + 1) * double(j + 2));
                    p(i, j) = p(j, i) = v;
                }
            }
            compare(basis, p, argv[arg], "dense_symmetric");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
