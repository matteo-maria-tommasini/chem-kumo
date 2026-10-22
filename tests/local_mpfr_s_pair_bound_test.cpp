#include "chem/kumo/local_mpfr_s_pair_bound.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace chem::kumo;
int main() {
    int checks=0;
    auto check=[&](bool ok){ ++checks; if(!ok) throw std::runtime_error("s-pair bound check failed at " + std::to_string(checks)); };
    GaussianPairIntermediate a{},b{};
    a.exponent=1.; b.exponent=1.; a.prefactor=1.; b.prefactor=1.;
    for(int d=0;d<3;++d) { a.polynomial[d][0]=1.; b.polynomial[d][0]=1.; }
    const double coincident=certified::s_pair_overlap_upper(a,b);
    check(coincident >= std::pow(std::acos(-1.)/2.,1.5));
    check(std::isfinite(coincident));
    b.center[0]=10.;
    const double separated=certified::s_pair_overlap_upper(a,b);
    check(separated>0. && separated<coincident);
    b.center[0]=100.;
    const double tiny=certified::s_pair_overlap_upper(a,b);
    check(tiny>0. && tiny<=std::numeric_limits<double>::denorm_min());
    b.center[0]=1.e150;
    check(certified::s_pair_overlap_upper(a,b)>0.);
    b.center[0]=0.; b.degree[0]=1;
    bool rejected=false;
    try { (void)certified::s_pair_overlap_upper(a,b); } catch(const std::invalid_argument&) {rejected=true;}
    check(rejected);
    std::cout << "checks="<<checks<<" s_pair_mpfr=OK status=OK\n";
}
