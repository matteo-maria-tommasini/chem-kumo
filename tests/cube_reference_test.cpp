#include <chem/kumo.hpp>
#include <Eigen/LU>
#include <array>
#include <limits>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace chem::kumo;
// Gaussian Cube headers round grid vectors to six decimal places. For the
// standard -3 grid, the true step is 1/6 bohr, not the printed 0.166667.
// Correct only components demonstrably within header rounding of 1/6.
// This is a reference-test compatibility adjustment, not a general Cube rule.
void recover_gaussian_standard_step(Grid& g) {
  constexpr double h = 1.0 / 6.0;
  constexpr double tolerance = 0.5e-6 + 1e-12;
  for (int a=0; a<3; ++a) {
    for (int c=0; c<3; ++c) {
      double& v=g.steps(c,a);
      if (c==a && std::abs(v-h)<tolerance) v=h;
    }
  }
}
int main(int argc, char** argv) {
 if(argc!=3 && argc!=4) {std::cerr<<"usage: kumo-cube-reference-test molecule.fchk reference.cube [total|alpha|beta|spin]\n"; return 2;}
 try {
  std::ifstream in(argv[2]); if(!in) throw std::runtime_error("cannot open reference cube");
  std::string line; std::getline(in,line); std::getline(in,line);
  int natoms, units; Grid g; if(!(in>>natoms>>g.origin[0]>>g.origin[1]>>g.origin[2]>>units)) throw std::runtime_error("invalid cube header");
  if(natoms<=0 || units!=1) throw std::runtime_error("unsupported cube header: expected positive atom count and one scalar field");
  for(int axis=0;axis<3;++axis) {int n; if(!(in>>n>>g.steps(0,axis)>>g.steps(1,axis)>>g.steps(2,axis)) || n<=0) throw std::runtime_error("invalid cube axis"); g.counts[axis]=n;}
  for(int a=0;a<natoms;++a) {int z; double charge,x,y,zpos; if(!(in>>z>>charge>>x>>y>>zpos)) throw std::runtime_error("invalid cube atom");}
  recover_gaussian_standard_step(g);
  auto mol=Molecule::from_fchk(argv[1]);
  const std::string component = argc==4 ? argv[3] : "total";
  Spin spin = Spin::Total;
  if (component=="alpha") spin=Spin::Alpha;
  else if (component=="beta") spin=Spin::Beta;
  else if (component=="spin") spin=Spin::Difference;
  else if (component!="total") throw std::runtime_error("invalid spin component: " + component);
  auto calculated=mol.densities_blocked(g,{spin},128).at(spin);
  const auto vals=calculated.values(); double max_abs=0, sum_sq=0, sum_abs=0, sum_ref=0, sum_val=0, sum_reference=0; std::size_t worst=0;
  double worst_ref=0, worst_val=0, worst_distance=0;
  std::array<double,3> shell_abs{}; std::array<double,3> shell_ref{};
  std::array<std::size_t,3> shell_count{};
  const int ny=g.counts[1], nz=g.counts[2];
  for(std::size_t i=0;i<vals.size();++i) {
    double ref;
    if(!(in>>ref) || !std::isfinite(ref)) throw std::runtime_error("invalid or truncated cube values");
    const int k=static_cast<int>(i%nz), j=static_cast<int>((i/nz)%ny), x=static_cast<int>(i/(static_cast<std::size_t>(ny)*nz));
    const Eigen::Vector3d r=g.origin+g.steps*Eigen::Vector3d(x,j,k);
    double nearest=std::numeric_limits<double>::infinity();
    for(const auto& atom:mol.atoms()) nearest=std::min(nearest,(r-atom.position).norm());
    const int shell=nearest<0.5?0:(nearest<2.0?1:2);
    const double d=vals[i]-ref;
    if(std::abs(d)>max_abs) {max_abs=std::abs(d);worst=i;worst_ref=ref;worst_val=vals[i];worst_distance=nearest;}
    sum_sq+=d*d;sum_abs+=std::abs(d);sum_ref+=std::abs(ref);sum_val+=vals[i];sum_reference+=ref;
    shell_abs[shell]+=std::abs(d);shell_ref[shell]+=std::abs(ref);++shell_count[shell];
  }
  double extra; if(in>>extra) throw std::runtime_error("extra cube values");
  const double rmse=std::sqrt(sum_sq/vals.size()), rel=sum_abs/sum_ref, volume=std::abs(g.steps.determinant());
  std::cout<<std::setprecision(12)<<"reference="<<argv[2]<<" points="<<vals.size()<<" max_abs="<<max_abs<<" rmse="<<rmse<<" relative_L1="<<rel<<" worst_index="<<worst<<" electron_sum_kumo="<<sum_val*volume<<" electron_sum_gaussian="<<sum_reference*volume<<'\n';
  const int wi=static_cast<int>(worst/(static_cast<std::size_t>(ny)*nz));
  const int wj=static_cast<int>((worst/nz)%ny), wk=static_cast<int>(worst%nz);
  const Eigen::Vector3d wr=g.origin+g.steps*Eigen::Vector3d(wi,wj,wk);
  std::cout<<"worst_point_ijk="<<wi<<","<<wj<<","<<wk<<" position_bohr="<<wr.transpose()
           <<" distance_nearest_nucleus_bohr="<<worst_distance<<" gaussian="<<worst_ref
           <<" kumo="<<worst_val<<" signed_difference="<<worst_val-worst_ref<<'\n';
  for(int s=0;s<3;++s) std::cout<<"distance_bin="<<(s==0?"[0,0.5)":s==1?"[0.5,2)":"[2,infinity)")
       <<" points="<<shell_count[s]<<" relative_L1="<<(shell_ref[s]>0?shell_abs[s]/shell_ref[s]:0)<<'\n';
  // A gross electron-count mismatch almost always indicates a valence-only
  // Gaussian cube instead of FDensity=SCF.
  if (mol.atoms().empty()) throw std::runtime_error("empty molecule");
  if (std::abs(sum_reference*volume-sum_val*volume)>0.5) {
    std::cerr<<"FAIL: reference electron count differs substantially from KUMO; check reference Cube component and full-density options\n";
    return 1;
  }
  if(max_abs>2e-3 || rel>1e-5) {std::cerr<<"FAIL: reference tolerance exceeded\n";return 1;}
  return 0;
 } catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
