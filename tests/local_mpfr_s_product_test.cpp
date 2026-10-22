#include "chem/kumo/local_mpfr_s_product.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace chem::kumo::certified;
int main() {
    int checks=0;
    auto check=[&](bool x){++checks;if(!x) throw std::runtime_error("s product test check "+std::to_string(checks));};
    SPrimitive a{1.,2.,{0.,0.,0.}},b{1.,3.,{0.,0.,0.}};
    auto x=s_primitive_product_intervals(a,b);
    check(x.exponent[0]<=2. && x.exponent[1]>=2.);
    check(x.reduced_exponent[0]<=0.5 && x.reduced_exponent[1]>=0.5);
    check(x.amplitude[0]<=6. && x.amplitude[1]>=6.);
    for(auto c:x.center) check(c[0]<=0. && c[1]>=0.);
    b.center[0]=10.; x=s_primitive_product_intervals(a,b);
    check(x.center[0][0]<=5. && x.center[0][1]>=5.);
    const double reference=6.*std::exp(-50.);
    check(x.amplitude[0]<=reference*(1.+1.e-14) && x.amplitude[1]>=reference*(1.-1.e-14));
    b.center[0]=100.; x=s_primitive_product_intervals(a,b);
    check(x.amplitude[1]>0. && x.amplitude[1]<=std::numeric_limits<double>::denorm_min());
    b.center[0]=1.e150; x=s_primitive_product_intervals(a,b);
    check(x.amplitude[1]>0.);
    b.center[0]=0.; b.exponent=std::numeric_limits<double>::denorm_min();
    x=s_primitive_product_intervals(a,b);
    check(x.exponent[0]>0. && x.amplitude[1]>0.);
    b.exponent=0.;
    bool invalid=false;
    try {(void)s_primitive_product_intervals(a,b);}catch(const std::invalid_argument&){invalid=true;}
    check(invalid);
    std::cout<<"checks="<<checks<<" s_product_mpfr=OK status=OK\n";
}
