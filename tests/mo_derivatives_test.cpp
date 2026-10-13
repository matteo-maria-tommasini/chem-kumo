#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    constexpr double h = 1.e-4;
    for (int file = 1; file < argc; ++file) {
        const auto mol = chem::kumo::Molecule::from_fchk(argv[file]);
        double max_gradient = 0.0, max_hessian = 0.0;
        for (const auto channel : {chem::kumo::Spin::Alpha, chem::kumo::Spin::Beta}) {
            const auto count = mol.orbital_count(channel);
            if (count == 0) throw std::runtime_error("missing orbitals");
            const std::size_t indices[] = {1, std::min<std::size_t>(3, count),
                std::min<std::size_t>(static_cast<std::size_t>(mol.occupied_count(channel)) + 1, count)};
            for (const auto index : indices) {
                for (int t = 0; t < 8; ++t) {
                    const Eigen::Vector3d r(0.41 + 0.33*t, -0.37 + 0.17*t, 0.23 - 0.13*t);
                    const auto exact = mol.orbital_derivatives(r, index, channel);
                    if (!std::isfinite(exact.value) || !exact.gradient.allFinite() ||
                        !exact.hessian.allFinite() ||
                        (exact.hessian - exact.hessian.transpose()).cwiseAbs().maxCoeff() > 1.e-13)
                        throw std::runtime_error("invalid analytic MO derivatives");
                    for (int a = 0; a < 3; ++a) {
                        Eigen::Vector3d rp = r, rm = r;
                        rp[a] += h; rm[a] -= h;
                        const auto plus = mol.orbital_derivatives(rp, index, channel);
                        const auto minus = mol.orbital_derivatives(rm, index, channel);
                        const double fd1 = (plus.value - minus.value)/(2*h);
                        max_gradient = std::max(max_gradient, std::abs(fd1 - exact.gradient[a]));
                        for (int b = 0; b < 3; ++b) {
                            const double fd2 = (plus.gradient[b] - minus.gradient[b])/(2*h);
                            max_hessian = std::max(max_hessian, std::abs(fd2 - exact.hessian(b,a)));
                        }
                    }
                }
            }
        }
        std::cout << argv[file] << " grad=" << max_gradient
                  << " hessian=" << max_hessian << '\n';
        if (max_gradient > 1.e-3 || max_hessian > 1.e-2)
            throw std::runtime_error("MO spatial derivative regression failed");
        bool rejected = false;
        try { (void)mol.orbital_derivatives(Eigen::Vector3d::Zero(), 0); }
        catch (const std::out_of_range&) { rejected = true; }
        if (!rejected) throw std::runtime_error("zero MO index accepted");
        rejected = false;
        try { (void)mol.orbital_derivatives(Eigen::Vector3d::Zero(), 1, chem::kumo::Spin::Total); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("invalid MO channel accepted");
    }
}
