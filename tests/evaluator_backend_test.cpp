#include "chem/kumo.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace chem::kumo::detail {
struct AoMoTestAccess {
    static auto run(const Molecule& m,
        const Eigen::Ref<const Eigen::Matrix<double,3,Eigen::Dynamic>>& points,
        const std::vector<Spin>& components, std::size_t block, int order) {
        return m.density_derivatives_mo(points,components,block,order);
    }
};
}
using namespace chem::kumo;
namespace {
void close(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b) ||
        std::abs(a-b)>2e-11*std::max({1.0,std::abs(a),std::abs(b)}))
        throw std::runtime_error("AO_MO Evaluator disagrees with direct MO kernel");
}
void run(const Molecule& m) {
    Grid grid;
    grid.origin=Eigen::Vector3d(0.21,-0.16,0.34);
    grid.steps << 0.23,0.01,0.0, -0.02,0.19,0.03, 0.01,0.0,0.17;
    grid.counts=Eigen::Vector3i(3,2,3);
    Eigen::Matrix<double,3,Eigen::Dynamic> points(3,18);
    for (int i=0;i<3;++i) for(int j=0;j<2;++j) for(int k=0;k<3;++k)
        points.col((i*2+j)*3+k)=grid.origin+grid.steps*Eigen::Vector3d(i,j,k);
    const std::array<DensityDescriptor,4> ds{density::total,density::alpha,density::beta,density::spin};
    const std::array<Spin,4> spins{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference};
    for (std::size_t block: {std::size_t{1},std::size_t{7},std::size_t{32}}) {
        Evaluator mo(m,{.threads=1,.block_size=block,.backend=Backend::AO_MO});
        Evaluator ao(m,{.threads=1,.block_size=block,.backend=Backend::AO_P});
        for(int order: {0,1,2}) {
            Request req;
            for(auto d:ds) {
                req.add(d.value());
                if(order>=1) req.add(d.gradient()).add(d.tau());
                if(order>=2) req.add(d.hessian()).add(d.laplacian());
            }
            const auto result=mo.evaluate(grid,req);
            if (result.backend_used() != Backend::AO_MO)
                throw std::runtime_error("explicit AO_MO backend metadata mismatch");
            const auto automatic = Evaluator(m,{.threads=1,.block_size=block,.backend=Backend::AUTO}).evaluate(grid,req);
            if (automatic.backend_used() != Backend::AO_MO)
                throw std::runtime_error("AUTO did not select eligible AO_MO");
            const auto reference=detail::AoMoTestAccess::run(m,points,
                {spins.begin(),spins.end()},block,order);
            for(std::size_t c=0;c<4;++c) {
                const auto& ref=reference.at(spins[c]);
                for(std::size_t p=0;p<18;++p) {
                    close(result.get(ds[c].value()).component()[p],ref[p].value);
                    if(order>=1) {
                        for(int a=0;a<3;++a)
                            close(result.get(ds[c].gradient()).component(a)[p],ref[p].gradient[a]);
                        close(result.get(ds[c].tau()).component()[p],ref[p].tau);
                    }
                    if(order>=2) {
                        constexpr int rows[6]={0,0,0,1,1,2};
                        constexpr int cols[6]={0,1,2,1,2,2};
                        for(int a=0;a<6;++a)
                            close(result.get(ds[c].hessian()).component(a)[p],ref[p].hessian(rows[a],cols[a]));
                        close(result.get(ds[c].laplacian()).component()[p],ref[p].laplacian);
                    }
                }
            }
            if(order==2) {
                const auto ao_result=ao.evaluate(grid,req);
                if (ao_result.backend_used() != Backend::AO_P)
                    throw std::runtime_error("explicit AO_P backend metadata mismatch");
                if(!ao_result.contains(density::total.hessian()))
                    throw std::runtime_error("AO_P backend regression");
            }
        }
    }
    bool rejected=false;
    try { Evaluator invalid(m,{.threads=1,.backend=static_cast<Backend>(99)}); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected) throw std::runtime_error("invalid backend accepted");
}
}
int main(int argc,char** argv) {
    if(argc<2) return 2;
    try {
        for(int i=1;i<argc;++i) {
            auto m=Molecule::from_fchk(argv[i]);
            if(!m.electronic_structure_info().ao_mo_eligible())
                throw std::runtime_error("fixture not AO_MO eligible");
            run(m);
        }
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
