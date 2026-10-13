#include "chem/kumo.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
  if(argc!=3) return 2;
  using namespace chem::kumo;
  for(int n=1;n<argc;++n) {
    auto mol=Molecule::from_fchk(argv[n]);
    auto grid=Grid::around(mol.atoms(),2.0,0.65);
    double maxerr=0;
    for(auto size:{std::size_t(1),std::size_t(16),std::size_t(32),std::size_t(64),std::size_t(128),std::size_t(256)}) {
      auto ref=mol.densities(grid,{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference});
      auto blk=mol.densities_blocked(grid,{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference},size);
      for(auto component:{Spin::Total,Spin::Alpha,Spin::Beta,Spin::Difference}) {
        auto x=ref.at(component).values(),y=blk.at(component).values();
        for(std::size_t k=0;k<x.size();++k) maxerr=std::max(maxerr,std::abs(x[k]-y[k]));
      }
    }
    auto only=mol.densities_blocked(grid,{Spin::Total,Spin::Total},32);
    if(only.size()!=1) throw std::runtime_error("duplicate component not eliminated");
    bool threw=false;
    try { (void)mol.densities_blocked(grid,{Spin::Total},0); } catch(const std::invalid_argument&) { threw=true; }
    if(!threw || maxerr>1e-9) throw std::runtime_error("blocked evaluation regression");
    std::cout<<argv[n]<<" max_abs="<<maxerr<<"\n";
  }
}
