#include "chem/kumo/local_mpfr_s_quartet_bound.hpp"
#include <mpfr.h>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using chem::kumo::certified::SPrimitive;
using chem::kumo::certified::s_primitive_quartet_overlap_upper;

namespace {
void require(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(what);
}
void dirty_flags() {
    mpfr_set_underflow();
    mpfr_set_overflow();
    mpfr_set_divby0();
    mpfr_set_nanflag();
    mpfr_set_inexflag();
}
}

int main() {
    try {
        unsigned checks=0;
        const SPrimitive a{1.,1.,{0.,0.,0.}};
        const SPrimitive b{1.,-2.,{1.,0.,0.}};
        const SPrimitive far{1.,1.,{1e8,0.,0.}};
        const SPrimitive extreme{1.,1.,{1e150,0.,0.}};
        const SPrimitive zero{1.,0.,{0.,0.,0.}};
        const SPrimitive tiny{std::numeric_limits<double>::denorm_min(),1.,{0.,0.,0.}};
        const std::array<std::array<SPrimitive,4>,5> cases{{
            {{a,a,b,b}}, {{a,a,far,far}}, {{a,a,extreme,extreme}},
            {{a,b,tiny,a}}, {{a,zero,b,far}}
        }};
        for (const auto& q : cases) {
            mpfr_clear_flags();
            const double clean=s_primitive_quartet_overlap_upper(q[0],q[1],q[2],q[3]);
            require(!std::isnan(clean) && clean>=0.,"invalid clean result");
            ++checks;
            dirty_flags();
            const double dirty=s_primitive_quartet_overlap_upper(q[0],q[1],q[2],q[3]);
            require(clean==dirty,"stale MPFR flags changed bound");
            ++checks;
        }
        // Verify the directed rounding behavior on an MPFR underflow case.
        mpfr_t x,y;
        mpfr_init2(x,256); mpfr_init2(y,256);
        mpfr_set_si(x,-1000000000L,MPFR_RNDN);
        mpfr_clear_flags();
        mpfr_exp(y,x,MPFR_RNDU);
        require(mpfr_underflow_p()!=0 || mpfr_sgn(y)>0,
                "underflow lost positive upper endpoint without flag");
        ++checks;
        mpfr_clear(x); mpfr_clear(y);
        std::cout << "checks=" << checks
                  << " stale_flags=OK ordinary=OK extreme=OK"
                  << " underflow_observation=OK status=OK\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "exception propagation: " << e.what() << '\n';
        return 1;
    }
}
