#include "chem/kumo.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
int main(int argc,char** argv) {
  if(argc<2 || argc>4) {std::cerr<<"usage: kumo-block-benchmark file.fchk [spacing=0.3] [repeats=5]\n";return 2;}
  using namespace chem::kumo;
  auto mol=Molecule::from_fchk(argv[1]);
  const double spacing=argc>2?std::stod(argv[2]):0.3;
  const int repeats=argc>3?std::stoi(argv[3]):5;
  if(repeats<1 || repeats>100) return 2;
  auto grid=Grid::around(mol.atoms(),3.0,spacing);
  volatile double sink=0;
  auto measure=[&](auto fn){double best=1e100;for(int r=0;r<repeats;++r){auto t0=std::chrono::steady_clock::now();auto result=fn();sink=sink+result.at(Spin::Total).values()[0];auto t1=std::chrono::steady_clock::now();best=std::min(best,std::chrono::duration<double>(t1-t0).count());}return best;};
  auto ref=measure([&]{std::map<Spin,Field> result;
      result.emplace(Spin::Total, mol.density_pointwise(grid,Spin::Total));
      result.emplace(Spin::Difference, mol.density_pointwise(grid,Spin::Difference));
      return result;});
  std::cout<<"basis="<<mol.basis_size()<<" points="<<static_cast<long long>(grid.counts[0])*grid.counts[1]*grid.counts[2]<<" repeats="<<repeats<<"\n";
  std::cout<<"pointwise_s="<<ref<<"\n";
  for(auto b:{16,32,64,128,256}) {auto t=measure([&]{return mol.densities_blocked(grid,{Spin::Total,Spin::Difference},b);});std::cout<<"block_"<<b<<"_s="<<t<<" speedup="<<ref/t<<"\n";}
  std::cout<<"sink="<<sink<<"\n";
}
