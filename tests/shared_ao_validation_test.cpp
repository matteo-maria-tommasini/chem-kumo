#include "chem/kumo.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace chem::kumo;
namespace {
constexpr std::array<DensityComponent,4> components{
    DensityComponent::Total, DensityComponent::Alpha,
    DensityComponent::Beta, DensityComponent::Spin};
constexpr std::array<Spin,4> spins{
    Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference};
void close(double got, double expected, const char* name) {
    const double tolerance = 5e-11 * std::max({1.0, std::abs(got), std::abs(expected)});
    if (!std::isfinite(got) || !std::isfinite(expected) || std::abs(got-expected)>tolerance)
        throw std::runtime_error(std::string(name)+": shared AO mismatch: " +
                                 std::to_string(got)+" vs "+std::to_string(expected));
}
Grid make_grid() {
    Grid grid;
    grid.origin = Eigen::Vector3d(0.13,-0.42,0.31);
    grid.steps << 0.19, 0.02, 0.01,
                 -0.03, 0.17, 0.02,
                  0.01,-0.01, 0.23;
    grid.counts = Eigen::Vector3i(5,4,7); // 140 points, partial block for B=128
    return grid;
}
void validate(const Molecule& molecule, std::size_t block_size) {
    const Grid grid = make_grid();
    const std::size_t n = 140;
    Eigen::Matrix<double,3,Eigen::Dynamic> points(3,static_cast<Eigen::Index>(n));
    for (int i=0;i<grid.counts[0];++i)
        for (int j=0;j<grid.counts[1];++j)
            for (int k=0;k<grid.counts[2];++k) {
                const auto p=static_cast<std::size_t>((i*grid.counts[1]+j)*grid.counts[2]+k);
                points.col(static_cast<Eigen::Index>(p))=grid.origin+grid.steps*Eigen::Vector3d(i,j,k);
            }
    Request request;
    for (auto c:components) {
        const DensityDescriptor d{c};
        request.add(d.gradient()).add(d.hessian()).add(d.laplacian()).add(d.tau());
    }
    Evaluator evaluator(molecule,{.threads=1,.block_size=block_size, .backend=Backend::AO_P});
    const auto fused=evaluator.evaluate(grid,request);
    for (std::size_t c=0;c<components.size();++c) {
        const auto reference=molecule.density_derivatives_blocked(points,spins[c],block_size);
        const DensityDescriptor d{components[c]};
        const auto& grad=fused.get(d.gradient());
        const auto& hess=fused.get(d.hessian());
        const auto& lap=fused.get(d.laplacian());
        const auto& tau=fused.get(d.tau());
        constexpr std::array<int,6> rows{0,0,0,1,1,2};
        constexpr std::array<int,6> cols{0,1,2,1,2,2};
        for (std::size_t p=0;p<n;++p) {
            for (int a=0;a<3;++a) close(grad.component(a)[p],reference[p].gradient[a],"gradient");
            for (int a=0;a<6;++a) close(hess.component(a)[p],reference[p].hessian(rows[a],cols[a]),"hessian");
            close(lap.component()[p],reference[p].laplacian,"laplacian");
            close(tau.component()[p],reference[p].tau,"tau");
            close(lap.component()[p],hess.component(0)[p]+hess.component(3)[p]+hess.component(5)[p],"trace");
        }
    }
    // Verify spin algebra for every primitive, including the kinetic density.
    for (auto kind:{PropertyKind::Gradient,PropertyKind::Hessian,PropertyKind::Laplacian,PropertyKind::Tau}) {
        const auto& total=fused.get(Property{kind,DensityComponent::Total});
        const auto& alpha=fused.get(Property{kind,DensityComponent::Alpha});
        const auto& beta=fused.get(Property{kind,DensityComponent::Beta});
        const auto& spin=fused.get(Property{kind,DensityComponent::Spin});
        for (std::size_t a=0;a<total.components.size();++a)
            for (std::size_t p=0;p<n;++p) {
                close(total.component(a)[p],alpha.component(a)[p]+beta.component(a)[p],"total spin identity");
                close(spin.component(a)[p],alpha.component(a)[p]-beta.component(a)[p],"difference spin identity");
            }
    }
}
}
int main(int argc,char** argv) {
    if (argc!=3) return 2;
    try {
        for (int i=1;i<argc;++i) {
            const auto molecule=Molecule::from_fchk(argv[i]);
            for (std::size_t b:{std::size_t{1},std::size_t{7},std::size_t{128},std::size_t{256}})
                validate(molecule,b);
        }
    } catch(const std::exception& e) {
        std::cerr << "shared AO validation: " << e.what() << '\n';
        return 1;
    }
}
