#include "chem/kumo/local_mpfr_s_quartet_bound.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace chem::kumo::certified;
int main() {
  int checks=0;
  auto check=[&](bool ok){++checks;if(!ok)throw std::runtime_error("quartet check "+std::to_string(checks));};
  SPrimitive a{1.,1.,{0.,0.,0.}},b=a,c=a,d=a;
  const double pi=std::acos(-1.0);
  const double expected=std::pow(pi/4.0,1.5);
  check(s_primitive_quartet_overlap_upper(a,b,c,d)>=expected);
  check(s_primitive_quartet_overlap_upper(a,b,c,d)<expected*(1+1.e-12));
  d.coefficient=-2.;
  check(s_primitive_quartet_overlap_upper(a,b,c,d)>=2*expected);
  d.coefficient=0.;
  check(s_primitive_quartet_overlap_upper(a,b,c,d)==0.);
  d.coefficient=1.; d.center[0]=10.;
  check(s_primitive_quartet_overlap_upper(a,b,c,d)>0.);
  d.center[0]=1.e150;
  check(s_primitive_quartet_overlap_upper(a,b,c,d)>0.);
  d.exponent=0.;
  bool threw=false;
  try {(void)s_primitive_quartet_overlap_upper(a,b,c,d);}catch(const std::invalid_argument&){threw=true;}
  check(threw);
  std::cout<<"checks="<<checks<<" s_quartet_mpfr=OK status=OK\n";
}
