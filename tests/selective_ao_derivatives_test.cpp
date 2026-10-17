#include "chem/kumo.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace chem::kumo;
namespace {
void check(double x, double y) {
    if (!std::isfinite(x) || !std::isfinite(y) ||
        std::abs(x-y) > 5e-11 * std::max({1.0,std::abs(x),std::abs(y)}))
        throw std::runtime_error("selective AO derivative mismatch");
}
void validate(const char* path) {
    const auto molecule = Molecule::from_fchk(path);
    Grid grid;
    grid.origin = Eigen::Vector3d(0.17,-0.31,0.26);
    grid.steps << 0.16,0.01,0.0, -0.02,0.21,0.01, 0.01,-0.02,0.19;
    grid.counts = Eigen::Vector3i(3,4,5);
    for (std::size_t block : {std::size_t{1},std::size_t{7},std::size_t{128}}) {
        Evaluator evaluator(molecule, EvaluationOptions{1,block});
        Request first, second;
        for (auto component : {DensityComponent::Total,DensityComponent::Alpha,
                               DensityComponent::Beta,DensityComponent::Spin}) {
            DensityDescriptor d{component};
            first.add(d.gradient()).add(d.tau());
            second.add(d.gradient()).add(d.tau()).add(d.hessian());
        }
        const auto low = evaluator.evaluate(grid,first);
        const auto full = evaluator.evaluate(grid,second);
        for (auto component : {DensityComponent::Total,DensityComponent::Alpha,
                               DensityComponent::Beta,DensityComponent::Spin}) {
            DensityDescriptor d{component};
            for (auto key : {d.gradient(),d.tau()}) {
                const auto& a=low.get(key);
                const auto& b=full.get(key);
                for (std::size_t c=0;c<a.components.size();++c)
                    for (std::size_t i=0;i<a.component(c).size();++i)
                        check(a.component(c)[i],b.component(c)[i]);
            }
        }
    }
}
}
int main(int argc,char** argv) {
    if (argc != 3) return 2;
    try { validate(argv[1]); validate(argv[2]); }
    catch (const std::exception& e) {
        std::cerr << "selective AO derivatives: " << e.what() << '\n';
        return 1;
    }
}
