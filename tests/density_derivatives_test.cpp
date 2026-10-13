#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double a, double b, double absolute, double relative, const char* message) {
    require(std::isfinite(a) && std::isfinite(b) &&
            std::abs(a-b) <= absolute + relative*std::max(std::abs(a),std::abs(b)), message);
}
}
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    constexpr double step = 1.e-4;
    for (int file = 1; file < argc; ++file) {
        const auto molecule = chem::kumo::Molecule::from_fchk(argv[file]);
        double max_gradient = 0.0, max_hessian = 0.0;
        for (int t = 0; t < 9; ++t) {
            const Eigen::Vector3d r(0.41 + 0.29*t, -0.37 + 0.13*t, 0.23 - 0.17*t);
            const auto total = molecule.density_derivatives(r);
            const auto alpha = molecule.density_derivatives(r, chem::kumo::Spin::Alpha);
            const auto beta = molecule.density_derivatives(r, chem::kumo::Spin::Beta);
            const auto spin = molecule.density_derivatives(r, chem::kumo::Spin::Difference);
            for (const auto* f : {&total, &alpha, &beta, &spin}) {
                require(std::isfinite(f->value) && f->gradient.allFinite() &&
                        f->hessian.allFinite() && std::isfinite(f->tau), "nonfinite derivative");
                close(f->laplacian, f->hessian.trace(), 1.e-12, 1.e-12, "Laplacian mismatch");
                require((f->hessian-f->hessian.transpose()).cwiseAbs().maxCoeff()<1.e-12,
                        "asymmetric Hessian");
            }
            for (int a=0; a<3; ++a) {
                Eigen::Vector3d rp=r, rm=r;
                rp[a]+=step; rm[a]-=step;
                for (auto component : {chem::kumo::Spin::Total, chem::kumo::Spin::Alpha,
                                       chem::kumo::Spin::Beta, chem::kumo::Spin::Difference}) {
                    const auto exact=molecule.density_derivatives(r,component);
                    const auto plus=molecule.density_derivatives(rp,component);
                    const auto minus=molecule.density_derivatives(rm,component);
                    const double fd=(plus.value-minus.value)/(2*step);
                    max_gradient=std::max(max_gradient,std::abs(fd-exact.gradient[a]));
                    close(fd,exact.gradient[a],2.e-4,2.e-5,"density gradient finite difference");
                    for (int b=0; b<3; ++b) {
                        const double second=(plus.gradient[b]-minus.gradient[b])/(2*step);
                        max_hessian=std::max(max_hessian,std::abs(second-exact.hessian(b,a)));
                        close(second,exact.hessian(b,a),3.e-3,3.e-5,"density Hessian finite difference");
                    }
                }
            }
            // All fields, including tau, must obey the same spin decomposition.
            close(total.value,alpha.value+beta.value,1.e-10,1.e-12,"total density split");
            close(spin.value,alpha.value-beta.value,1.e-10,1.e-12,"spin density split");
            close(total.laplacian,alpha.laplacian+beta.laplacian,1.e-9,1.e-12,"Laplacian split");
            close(spin.laplacian,alpha.laplacian-beta.laplacian,1.e-9,1.e-12,"spin Laplacian split");
            close(total.tau,alpha.tau+beta.tau,1.e-10,1.e-12,"tau split");
            close(spin.tau,alpha.tau-beta.tau,1.e-10,1.e-12,"spin tau split");
            require((total.gradient-alpha.gradient-beta.gradient).norm()<1.e-9,"gradient split");
            require((spin.gradient-alpha.gradient+beta.gradient).norm()<1.e-9,"spin gradient split");
            require((total.hessian-alpha.hessian-beta.hessian).norm()<1.e-8,"Hessian split");
            require((spin.hessian-alpha.hessian+beta.hessian).norm()<1.e-8,"spin Hessian split");
            // Independent pointwise density implementation, not derivative code.
            chem::kumo::Grid one;
            one.origin=r;
            one.counts=Eigen::Vector3i::Ones();
            close(total.value,molecule.density_pointwise(one)(0,0,0),1.e-10,1.e-12,
                  "density value differs from reference path");
        }
        bool rejected=false;
        try { (void)molecule.density_derivatives(Eigen::Vector3d::Constant(NAN)); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"nonfinite point accepted");
        std::cout << argv[file] << " grad=" << max_gradient << " hessian=" << max_hessian << '\n';
    }
}
