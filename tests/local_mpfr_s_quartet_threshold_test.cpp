#include "chem/kumo/local_mpfr_s_quartet_bound.hpp"
#include <mpfr.h>
#include <cmath>
#include <initializer_list>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>

using chem::kumo::certified::SPrimitive;
using chem::kumo::certified::s_primitive_quartet_overlap_upper;

static void require(bool yes, const char* why) {
    if (!yes) throw std::runtime_error(why);
}

// Independent four-center identity, evaluated as a 1024-bit interval.
// In particular, the reference does not construct Gaussian product centers.
static void reference_bounds(const SPrimitive (&s)[4], mpfr_t lower, mpfr_t upper) {
    constexpr mpfr_prec_t bits=1024;
    mpfr_t alpha[4], sum, decay_lo, decay_hi, delta, sq, pair_lo, pair_hi;
    mpfr_t tmp, pi, ratio, volume, coeff, result;
    for(auto& a:alpha) mpfr_init2(a,bits);
    mpfr_inits2(bits,sum,decay_lo,decay_hi,delta,sq,pair_lo,pair_hi,
               tmp,pi,ratio,volume,coeff,result,(mpfr_ptr)0);
    mpfr_set_zero(sum,0);
    for(int i=0;i<4;++i) {
        mpfr_set_d(alpha[i],s[i].exponent,MPFR_RNDN);
        mpfr_add(sum,sum,alpha[i],MPFR_RNDU);
    }
    // All distances here are exact binary64 differences for the chosen
    // boundary inputs (zero and one variable coordinate); calculate
    // independent lower and upper enclosures nevertheless.
    mpfr_set_zero(decay_lo,0); mpfr_set_zero(decay_hi,0);
    for(int i=0;i<4;++i) for(int j=i+1;j<4;++j) {
        mpfr_set_zero(pair_lo,0); mpfr_set_zero(pair_hi,0);
        for(int k=0;k<3;++k) {
            mpfr_set_d(delta,s[i].center[k],MPFR_RNDN);
            mpfr_sub_d(delta,delta,s[j].center[k],MPFR_RNDN);
            mpfr_sqr(sq,delta,MPFR_RNDU);
            mpfr_add(pair_hi,pair_hi,sq,MPFR_RNDU);
            mpfr_sqr(sq,delta,MPFR_RNDD);
            mpfr_add(pair_lo,pair_lo,sq,MPFR_RNDD);
        }
        mpfr_mul(pair_lo,pair_lo,alpha[i],MPFR_RNDD);
        mpfr_mul(pair_lo,pair_lo,alpha[j],MPFR_RNDD);
        mpfr_mul(pair_hi,pair_hi,alpha[i],MPFR_RNDU);
        mpfr_mul(pair_hi,pair_hi,alpha[j],MPFR_RNDU);
        mpfr_add(decay_lo,decay_lo,pair_lo,MPFR_RNDD);
        mpfr_add(decay_hi,decay_hi,pair_hi,MPFR_RNDU);
    }
    mpfr_div(decay_lo,decay_lo,sum,MPFR_RNDD);
    mpfr_div(decay_hi,decay_hi,sum,MPFR_RNDU);
    mpfr_const_pi(pi,MPFR_RNDU);
    mpfr_div(ratio,pi,sum,MPFR_RNDU);
    mpfr_sqrt(volume,ratio,MPFR_RNDU);
    mpfr_mul(volume,volume,ratio,MPFR_RNDU);
    mpfr_set_ui(coeff,1,MPFR_RNDN);
    for(const auto& x:s) mpfr_mul_d(coeff,coeff,std::abs(x.coefficient),MPFR_RNDU);
    mpfr_mul(volume,volume,coeff,MPFR_RNDU);
    mpfr_neg(tmp,decay_lo,MPFR_RNDU);
    mpfr_exp(result,tmp,MPFR_RNDU);
    mpfr_mul(upper,volume,result,MPFR_RNDU);
    // A deliberately weak but valid lower bound: zero. The useful
    // comparison is to the independently rounded-up reference.
    mpfr_set_zero(lower,0);
    for(auto& a:alpha) mpfr_clear(a);
    mpfr_clears(sum,decay_lo,decay_hi,delta,sq,pair_lo,pair_hi,
                tmp,pi,ratio,volume,coeff,result,(mpfr_ptr)0);
}

int main() {
    try {
        mpfr_t lo,hi; mpfr_init2(lo,1024); mpfr_init2(hi,1024);
        const double shifts[]={0,1,10,20,30,35,38,39,40,41,50,100,1e6};
        std::size_t checked=0, subnormal=0, infinite=0;
        for(double x:shifts) for(double scale:{0.5,1.0,2.0}) {
            SPrimitive s[4]={
                {scale,1.0,{0,0,0}}, {scale,-1.0,{0,0,0}},
                {scale,1.0,{x,0,0}}, {scale,1.0,{x,0,0}}
            };
            const double bound=s_primitive_quartet_overlap_upper(s[0],s[1],s[2],s[3]);
            reference_bounds(s,lo,hi);
            const double ref=mpfr_get_d(hi,MPFR_RNDU);
            require(!std::isnan(bound) && bound>=0,"invalid bound");
            // Reference upper is independently conservative; comparing
            // against it is intentionally stronger than necessary.
            require(bound>=ref,"bound smaller than independent upper reference");
            if(std::fpclassify(bound)==FP_SUBNORMAL) ++subnormal;
            if(std::isinf(bound)) ++infinite;
            ++checked;
        }
        mpfr_clear(lo); mpfr_clear(hi);
        std::cout<<"cases="<<checked<<" subnormal="<<subnormal
                 <<" infinite="<<infinite<<" failures=0 status=OK\n";
    } catch(const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1;
    }
}
