#include "chem/kumo/local_mpfr_s_quartet_bound.hpp"
#include <mpfr.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
using chem::kumo::certified::SPrimitive;
using chem::kumo::certified::s_primitive_quartet_overlap_upper;
namespace {
struct M {
    mpfr_t x;
    M() { mpfr_init2(x,1024); }
    explicit M(double v):M() { mpfr_set_d(x,v,MPFR_RNDN); }
    ~M(){mpfr_clear(x);}
    M(const M&)=delete; M& operator=(const M&)=delete;
};
// Independent positive-term identity: sum_{i<j} alpha_i alpha_j |Ri-Rj|^2 / sum_i alpha_i.
// All operations below round upwards except the final negative exponent.
double reference_upper(const std::array<SPrimitive,4>& p) {
    for(const auto& v:p) if(v.coefficient==0.) return 0.;
    M sum(0.),numerator(0.),tmp,delta,square,term,weight,decay,negative,exponential,pi,ratio,root,volume,answer;
    for(const auto& v:p) {
        M a(v.exponent);
        mpfr_add(sum.x,sum.x,a.x,MPFR_RNDD); // lower denominator
    }
    for(int i=0;i<4;++i) for(int j=i+1;j<4;++j) {
        M ai(p[i].exponent),aj(p[j].exponent),dist2(0.);
        for(int k=0;k<3;++k) {
            M ri(p[i].center[k]),rj(p[j].center[k]);
            // Difference magnitude bounded by outward interval.
            M lo,hi,mag;
            mpfr_sub(lo.x,ri.x,rj.x,MPFR_RNDD);
            mpfr_sub(hi.x,ri.x,rj.x,MPFR_RNDU);
            mpfr_abs(lo.x,lo.x,MPFR_RNDU);
            mpfr_abs(hi.x,hi.x,MPFR_RNDU);
            if(mpfr_cmp(lo.x,hi.x)>0) mpfr_set(mag.x,lo.x,MPFR_RNDU);
            else mpfr_set(mag.x,hi.x,MPFR_RNDU);
            mpfr_sqr(square.x,mag.x,MPFR_RNDU);
            mpfr_add(dist2.x,dist2.x,square.x,MPFR_RNDU);
        }
        mpfr_mul(weight.x,ai.x,aj.x,MPFR_RNDU);
        mpfr_mul(term.x,weight.x,dist2.x,MPFR_RNDU);
        mpfr_add(numerator.x,numerator.x,term.x,MPFR_RNDU);
    }
    // For an upper bound on exp(-decay), decay must be rounded DOWN.
    // Compute a separate lower numerator with directed arithmetic.
    M numerator_lo(0.), dist_lo, diff_lo, sq_lo, w_lo, t_lo;
    for(int i=0;i<4;++i) for(int j=i+1;j<4;++j) {
        M ai(p[i].exponent),aj(p[j].exponent),distance(0.);
        for(int k=0;k<3;++k) {
            M ri(p[i].center[k]),rj(p[j].center[k]);
            mpfr_sub(diff_lo.x,ri.x,rj.x,MPFR_RNDN); // exact for binary64 coordinates at 1024 bits
            mpfr_sqr(sq_lo.x,diff_lo.x,MPFR_RNDD);
            mpfr_add(distance.x,distance.x,sq_lo.x,MPFR_RNDD);
        }
        mpfr_mul(w_lo.x,ai.x,aj.x,MPFR_RNDD);
        mpfr_mul(t_lo.x,w_lo.x,distance.x,MPFR_RNDD);
        mpfr_add(numerator_lo.x,numerator_lo.x,t_lo.x,MPFR_RNDD);
    }
    M sum_hi(0.);
    for(const auto& v:p){M a(v.exponent);mpfr_add(sum_hi.x,sum_hi.x,a.x,MPFR_RNDU);}
    mpfr_div(decay.x,numerator_lo.x,sum_hi.x,MPFR_RNDD);
    mpfr_neg(negative.x,decay.x,MPFR_RNDU);
    mpfr_exp(exponential.x,negative.x,MPFR_RNDU);
    mpfr_const_pi(pi.x,MPFR_RNDU);
    mpfr_div(ratio.x,pi.x,sum.x,MPFR_RNDU);
    mpfr_sqrt(root.x,ratio.x,MPFR_RNDU);
    mpfr_mul(volume.x,ratio.x,root.x,MPFR_RNDU);
    mpfr_mul(answer.x,volume.x,exponential.x,MPFR_RNDU);
    for(const auto& v:p) mpfr_mul_d(answer.x,answer.x,std::abs(v.coefficient),MPFR_RNDU);
    return mpfr_get_d(answer.x,MPFR_RNDU);
}
}
int main() {
    std::mt19937_64 rng(0x16d6c5ULL);
    std::uniform_real_distribution<double> uniform(-1.,1.);
    std::uniform_int_distribution<int> expo(-250,250);
    std::size_t cases=0,failed=0,finite=0,subnormal=0,infinite=0;
    auto check=[&](const std::array<SPrimitive,4>& p) {
        ++cases;
        double got=s_primitive_quartet_overlap_upper(p[0],p[1],p[2],p[3]);
        double ref=reference_upper(p);
        if(std::isfinite(got)) ++finite; else ++infinite;
        if(got>0 && got<std::numeric_limits<double>::min()) ++subnormal;
        if(!(got>=ref) || std::isnan(got)) {
            ++failed;
            if(failed<=8) std::cerr<<"FAIL case="<<cases<<" got="<<got<<" ref_upper="<<ref<<'\n';
        }
    };
    for(int n=0;n<3000;++n) {
        std::array<SPrimitive,4> p;
        for(auto& v:p) {
            v.exponent=std::ldexp(1.+0.9*std::abs(uniform(rng)),expo(rng));
            v.coefficient=(uniform(rng)<0?-1.:1.)*std::ldexp(1.+std::abs(uniform(rng)),expo(rng)/4);
            for(double& x:v.center) x=uniform(rng)*std::ldexp(1.,expo(rng)/3);
        }
        if(n%7==0) for(int j=1;j<4;++j) p[j].center=p[0].center;
        if(n%11==0) p[3].center[0]=1.e150;
        if(n%13==0) p[2].exponent=std::numeric_limits<double>::denorm_min();
        check(p);
    }
    std::cout<<"cases="<<cases<<" failures="<<failed<<" finite="<<finite
             <<" subnormal="<<subnormal<<" infinite="<<infinite
             <<" reference_precision=1024 status="<<(failed?"FAIL":"OK")<<'\n';
    return failed?1:0;
}
