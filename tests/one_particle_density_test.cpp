#include "chem/kumo/density.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;

static double maxabs(const Eigen::MatrixXd& m) {
    return m.cwiseAbs().maxCoeff();
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    try {
        for (int arg=1; arg<argc; ++arg) {
            const auto m = Molecule::from_fchk(argv[arg]);
            FormattedCheckpoint f(argv[arg]);
            const Eigen::MatrixXd S = BasisEvaluator{BasisSet(f)}.overlap_matrix();
            const auto pt = m.density_matrix(Spin::Total);
            const auto ps = m.density_matrix(Spin::Difference);
            const auto pa = m.density_matrix(Spin::Alpha);
            const auto pb = m.density_matrix(Spin::Beta);
            if (maxabs(pa + pb - pt) > 2e-14 ||
                maxabs(pa - pb - ps) > 2e-14)
                throw std::runtime_error("spin decomposition mismatch");
            const auto verify = [&](const Eigen::MatrixXd& p, double electrons, const char* label) {
                if (!p.allFinite() || maxabs(p-p.transpose()) > 1e-12)
                    throw std::runtime_error("nonfinite or asymmetric 1-RDM");
                const double trace = (p.cwiseProduct(S)).sum();
                const double idempotency = maxabs(p*S*p-p);
                // The SCF density is read from FCHK; it is not overwritten
                // by an occupied-orbital reconstruction. Tolerances account
                // for printed FCHK precision and AO metric conditioning.
                if (std::abs(trace-electrons) > 2e-4)
                    throw std::runtime_error("electron count mismatch");
                if (m.electronic_structure_info().ao_mo_eligible() && idempotency > 2e-4)
                    throw std::runtime_error("metric idempotency mismatch");
                std::cout << argv[arg] << " spin=" << label
                          << " electrons=" << std::setprecision(12) << trace
                          << " idempotency_max=" << idempotency << '\n';
            };
            verify(pa, m.occupied_count(Spin::Alpha), "alpha");
            verify(pb, m.occupied_count(Spin::Beta), "beta");
            const Eigen::Vector3d point(0.29,-0.41,0.37);
            const Eigen::VectorXd ao = BasisEvaluator{BasisSet(f)}.evaluate(point);
            for (Spin spin : {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference}) {
                const Eigen::MatrixXd p = m.density_matrix(spin);
                const double direct = ao.dot(p*ao);
                const double reference = m.density_derivatives(point,spin).value;
                if (std::abs(direct-reference) > 1e-9*std::max({1.0,std::abs(direct),std::abs(reference)}))
                    throw std::runtime_error("real-space density mismatch");
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
