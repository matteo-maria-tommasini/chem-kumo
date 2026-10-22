#include "chem/kumo/local_mpfr_s_quartet_bound.hpp"
#include <mpfr.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using chem::kumo::certified::SPrimitive;
using chem::kumo::certified::s_primitive_quartet_overlap_upper;

static void require(bool ok, const char* label) {
    if (!ok) throw std::runtime_error(label);
}

int main() {
    try {
        unsigned checks=0;
        // Directly demonstrate that MPFR's exponent range can underflow
        // even at 256-bit precision; no Gaussian-bound inference is made.
        mpfr_t x,y;
        mpfr_init2(x,256); mpfr_init2(y,256);
        mpfr_set_si(x,-1000000000L,MPFR_RNDN);
        mpfr_clear_flags();
        mpfr_exp(y,x,MPFR_RNDU);
        const bool underflow=mpfr_underflow_p()!=0;
        require(underflow || mpfr_sgn(y)>0,"MPFR exponential lost positive upper endpoint");
        ++checks;
        mpfr_clear_flags();
        mpfr_set_ui(x,1,MPFR_RNDN);
        mpfr_exp(y,x,MPFR_RNDU);
        require(mpfr_number_p(y) && !mpfr_underflow_p() && !mpfr_overflow_p(),"ordinary MPFR exponential flags");
        ++checks;
        mpfr_clear(x); mpfr_clear(y);

        const SPrimitive a{1.0,1.0,{0.0,0.0,0.0}};
        const SPrimitive near{1.0,1.0,{1.0,0.0,0.0}};
        const SPrimitive far{1.0,1.0,{1e8,0.0,0.0}};
        const SPrimitive extreme{1.0,1.0,{1e150,0.0,0.0}};
        const SPrimitive zero{1.0,0.0,{0.0,0.0,0.0}};
        const double ordinary=s_primitive_quartet_overlap_upper(a,a,near,near);
        require(std::isfinite(ordinary) && ordinary>0.0,"ordinary quartet"); ++checks;
        const double tiny=s_primitive_quartet_overlap_upper(a,a,far,far);
        require(tiny>=std::numeric_limits<double>::denorm_min() && !std::isnan(tiny),"far quartet lost conservative fallback"); ++checks;
        const double extreme_bound=s_primitive_quartet_overlap_upper(a,a,extreme,extreme);
        require(extreme_bound>=std::numeric_limits<double>::denorm_min() && !std::isnan(extreme_bound),"extreme quartet invalid"); ++checks;
        require(s_primitive_quartet_overlap_upper(a,a,zero,a)==0.0,"exact zero coefficient"); ++checks;
        std::cout << "checks=" << checks << " mpfr_underflow_observed=" << underflow
                  << " far_upper=" << tiny << " extreme_upper=" << extreme_bound
                  << " status=OK\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "MPFR exception safety: " << e.what() << '\n';
        return 1;
    }
}
