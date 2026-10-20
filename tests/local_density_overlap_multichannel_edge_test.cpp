#include "chem/kumo/local_density_overlap_multichannel.hpp"
#include "chem/kumo/local_density_overlap_symmetric_paired.hpp"
#include "chem/kumo/fchk.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Matrix = Eigen::MatrixXd;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void compare(const Matrix& actual, const Matrix& expected, const char* label) {
    require(actual.rows() == expected.rows() && actual.cols() == expected.cols(), "output shape mismatch");
    require(actual.allFinite(), "nonfinite output");
    const double scale = std::max(1.0, expected.cwiseAbs().maxCoeff());
    const double error = (actual - expected).cwiseAbs().maxCoeff();
    const double symmetry = (actual - actual.transpose()).cwiseAbs().maxCoeff();
    std::cout << " case=" << label << " error=" << error << " symmetry=" << symmetry;
    require(error <= 1e-11 * scale, "reference mismatch");
    require(symmetry <= 1e-12 * scale, "asymmetric output");
}

template<class F>
void must_throw(F&& f, const char* label) {
    bool threw = false;
    try { f(); } catch (const std::invalid_argument&) { threw = true; }
    if (!threw) throw std::runtime_error(std::string("expected invalid_argument: ") + label);
    std::cout << " rejected=" << label;
}
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        using namespace chem::kumo;
        std::cout << std::scientific << std::setprecision(12);
        for (int arg = 1; arg < argc; ++arg) {
            const BasisSet basis{FormattedCheckpoint(argv[arg])};
            const Eigen::Index n = static_cast<Eigen::Index>(basis.function_count());
            require(n >= 2, "basis too small");
            const Matrix zero = Matrix::Zero(n, n);
            Matrix dense(n, n);
            for (Eigen::Index i = 0; i < n; ++i)
                for (Eigen::Index j = 0; j < n; ++j)
                    dense(i, j) = 0.2 * std::cos(0.17 * static_cast<double>(i + j + 2))
                                + 0.1 * std::sin(0.13 * static_cast<double>((i + 1) * (j + 1)));
            dense = Matrix(0.5 * (dense + dense.transpose()).eval());
            std::cout << "system=" << argv[arg] << " nao=" << n;

            const auto empty = local_density_overlap_multichannel(basis, {});
            require(empty.matrices.empty(), "empty channel result is not empty");
            require(empty.statistics.evaluated_quartets == 0, "empty channels evaluated quartets");
            require(empty.statistics.active_channels == 0, "empty channel active count");
            std::cout << " empty=ok";

            const auto all_zero = local_density_overlap_multichannel(basis, {zero, zero});
            require(all_zero.matrices.size() == 2, "all-zero channel count mismatch");
            require(all_zero.statistics.active_channels == 0, "all-zero channels marked active");
            require(all_zero.statistics.evaluated_quartets == 0, "all-zero channels evaluated quartets");
            require(all_zero.statistics.prepared_pair_instances == 0, "all-zero channels prepared pairs");
            compare(all_zero.matrices[0], zero, "all_zero_first");
            compare(all_zero.matrices[1], zero, "all_zero_second");

            Matrix tiny = Matrix::Zero(n, n);
            tiny(0, 0) = 1e-200;
            const auto tiny_result = local_density_overlap_multichannel(basis, {tiny});
            require(tiny_result.statistics.active_channels == 1, "tiny nonzero channel marked inactive");
            require(tiny_result.statistics.evaluated_quartets > 0, "tiny nonzero channel skipped quartets");
            std::cout << " tiny_nonzero=active";

            const auto one = local_density_overlap_multichannel(basis, {dense});
            require(one.matrices.size() == 1, "single channel count mismatch");
            const Matrix reference = local_density_overlap_symmetric_paired(basis, dense);
            compare(one.matrices[0], reference, "dense_single");

            const auto multiple = local_density_overlap_multichannel(basis, {dense, dense, zero});
            require(multiple.matrices.size() == 3, "multiple channel count mismatch");
            require(multiple.statistics.evaluated_quartets == one.statistics.evaluated_quartets,
                    "quartets unexpectedly depend on channel count");
            require(multiple.statistics.active_channels == 2, "mixed channel active count mismatch");
            require(one.statistics.active_channels == 1, "single channel active count mismatch");
            compare(multiple.matrices[0], reference, "dense_first");
            compare(multiple.matrices[1], reference, "dense_duplicate");
            compare(multiple.matrices[2], zero, "zero");

            must_throw([&] { (void)local_density_overlap_multichannel(basis, {Matrix::Zero(n-1, n-1)}); }, "shape");
            Matrix asymmetric = dense;
            asymmetric(0, 1) += 0.01;
            must_throw([&] { (void)local_density_overlap_multichannel(basis, {asymmetric}); }, "asymmetry");
            Matrix nonfinite = dense;
            nonfinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
            must_throw([&] { (void)local_density_overlap_multichannel(basis, {nonfinite}); }, "nonfinite");
            std::cout << " quartets=" << one.statistics.evaluated_quartets << " OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "multichannel edge regression: " << e.what() << '\n';
        return 1;
    }
}
