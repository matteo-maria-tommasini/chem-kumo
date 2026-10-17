#include "chem/kumo.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace chem::kumo;
namespace {
volatile double sink = 0.0;
using Clock=std::chrono::steady_clock;
template<class F> double timed(F&& f, int repeats) {
    std::vector<double> seconds;
    for (int r=0;r<repeats;++r) {
        const auto begin=Clock::now();
        const double value=f();
        const auto end=Clock::now();
        sink = sink + value;
        seconds.push_back(std::chrono::duration<double>(end-begin).count());
    }
    std::sort(seconds.begin(),seconds.end());
    return seconds[seconds.size()/2];
}
}
int main(int argc,char** argv) {
    if (argc<2 || argc>5) {
        std::cerr << "usage: kumo-shared-ao-benchmark FILE.fchk[.gz] [points=512] [block=128] [repeats=5]\n";
        return 2;
    }
    try {
        const auto molecule=Molecule::from_fchk(argv[1]);
        const int count=argc>2?std::stoi(argv[2]):512;
        const int block=argc>3?std::stoi(argv[3]):128;
        const int repeats=argc>4?std::stoi(argv[4]):5;
        if(count<1 || block<1 || repeats<1) throw std::invalid_argument("positive parameters required");
        Grid grid;
        grid.origin=Eigen::Vector3d(-1.8,-1.2,-0.9);
        grid.steps=Eigen::Matrix3d::Identity()*0.12;
        grid.counts=Eigen::Vector3i(count,1,1);
        Eigen::Matrix<double,3,Eigen::Dynamic> points(3,count);
        for(int p=0;p<count;++p)
            points.col(p)=grid.origin+grid.steps*Eigen::Vector3d(p,0,0);
        Request request;
        for(auto c:{DensityComponent::Total,DensityComponent::Alpha,DensityComponent::Beta,DensityComponent::Spin}) {
            const DensityDescriptor d{c};
            request.add(d.gradient()).add(d.hessian()).add(d.laplacian()).add(d.tau());
        }
        Evaluator evaluator(molecule,{.threads=1,.block_size=static_cast<std::size_t>(block)});
        const auto old=[&]() {
            double checksum=0;
            for(auto s:{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference}) {
                const auto fields=molecule.density_derivatives_blocked(points,s,static_cast<std::size_t>(block));
                for(const auto& d:fields) checksum+=d.gradient.x()+d.tau+d.hessian(0,0);
            }
            return checksum;
        };
        const auto shared=[&]() {
            const auto fields=evaluator.evaluate(grid,request);
            double checksum=0;
            for(auto c:{DensityComponent::Total,DensityComponent::Alpha,DensityComponent::Beta,DensityComponent::Spin}) {
                const auto& grad=fields.get(Property{PropertyKind::Gradient,c});
                const auto& tau=fields.get(Property{PropertyKind::Tau,c});
                const auto& hess=fields.get(Property{PropertyKind::Hessian,c});
                for(int p=0;p<count;++p)
                    checksum+=grad.component(0)[p]+tau.component()[p]+hess.component(0)[p];
            }
            return checksum;
        };
        const double a=old(),b=shared();
        if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>5e-10*std::max({1.0,std::abs(a),std::abs(b)}))
            throw std::runtime_error("checksum mismatch");
        // Warm both implementations before timing. Timings are medians of repeated calls.
        (void)old(); (void)shared();
        const double old_time=timed(old,repeats);
        const double shared_time=timed(shared,repeats);
        std::cout << "points="<<count<<" block="<<block<<" repeats="<<repeats
                  <<" old_seconds="<<old_time<<" shared_seconds="<<shared_time
                  <<" speedup="<<old_time/shared_time<<" checksum="<<b<<'\n';
    } catch(const std::exception& e) {
        std::cerr << "shared AO benchmark: "<<e.what()<<'\n';
        return 1;
    }
}
