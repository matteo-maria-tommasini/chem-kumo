#pragma once
// Stage 03C.16d6c3. Directed MPFR enclosure of the s-s Gaussian product.
// All binary64 inputs denote their exact real values. Optional MPFR component.
#include <mpfr.h>
#include <array>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace chem::kumo::certified {
struct SPrimitive {
    double exponent;
    double coefficient;
    std::array<double,3> center;
};
struct SProductIntervals {
    std::array<double,2> exponent;
    std::array<double,2> reduced_exponent;
    std::array<std::array<double,2>,3> center;
    std::array<double,2> amplitude;
};
namespace sproduct_detail {
struct Num {
    mpfr_t v;
    Num() { mpfr_init2(v,256); }
    explicit Num(double d) : Num() { mpfr_set_d(v,d,MPFR_RNDN); }
    ~Num() { mpfr_clear(v); }
    Num(const Num&)=delete;
    Num& operator=(const Num&)=delete;
};
struct Range {
    Num lo,hi;
    Range()=default;
    explicit Range(double d) { mpfr_set_d(lo.v,d,MPFR_RNDN); mpfr_set_d(hi.v,d,MPFR_RNDN); }
    Range(const Range&)=delete;
    Range& operator=(const Range&)=delete;
};
inline void add(Range& z,const Range& x,const Range& y) {
    mpfr_add(z.lo.v,x.lo.v,y.lo.v,MPFR_RNDD);
    mpfr_add(z.hi.v,x.hi.v,y.hi.v,MPFR_RNDU);
}
inline void sub(Range& z,const Range& x,const Range& y) {
    mpfr_sub(z.lo.v,x.lo.v,y.hi.v,MPFR_RNDD);
    mpfr_sub(z.hi.v,x.hi.v,y.lo.v,MPFR_RNDU);
}
inline void mul(Range& z,const Range& x,const Range& y) {
    Num low,high,t;
    mpfr_mul(low.v,x.lo.v,y.lo.v,MPFR_RNDD);
    mpfr_mul(high.v,x.lo.v,y.lo.v,MPFR_RNDU);
    const mpfr_srcptr xs[2]={x.lo.v,x.hi.v};
    const mpfr_srcptr ys[2]={y.lo.v,y.hi.v};
    for(int i=0;i<2;++i) for(int j=0;j<2;++j) {
        mpfr_mul(t.v,xs[i],ys[j],MPFR_RNDD);
        if(mpfr_cmp(t.v,low.v)<0) mpfr_set(low.v,t.v,MPFR_RNDD);
        mpfr_mul(t.v,xs[i],ys[j],MPFR_RNDU);
        if(mpfr_cmp(t.v,high.v)>0) mpfr_set(high.v,t.v,MPFR_RNDU);
    }
    mpfr_set(z.lo.v,low.v,MPFR_RNDD);
    mpfr_set(z.hi.v,high.v,MPFR_RNDU);
}
inline void div_positive(Range& z,const Range& x,const Range& y) {
    if(mpfr_sgn(y.lo.v)<=0) throw std::domain_error("nonpositive interval divisor");
    Range reciprocal;
    Num one(1.0);
    mpfr_div(reciprocal.lo.v,one.v,y.hi.v,MPFR_RNDD);
    mpfr_div(reciprocal.hi.v,one.v,y.lo.v,MPFR_RNDU);
    mul(z,x,reciprocal);
}
inline void square(Range& z,const Range& x) {
    if(mpfr_sgn(x.lo.v)>=0) {
        mpfr_sqr(z.lo.v,x.lo.v,MPFR_RNDD);
        mpfr_sqr(z.hi.v,x.hi.v,MPFR_RNDU);
    } else if(mpfr_sgn(x.hi.v)<=0) {
        mpfr_sqr(z.lo.v,x.hi.v,MPFR_RNDD);
        mpfr_sqr(z.hi.v,x.lo.v,MPFR_RNDU);
    } else {
        mpfr_set_zero(z.lo.v,0);
        Num l,r;
        mpfr_sqr(l.v,x.lo.v,MPFR_RNDU);
        mpfr_sqr(r.v,x.hi.v,MPFR_RNDU);
        mpfr_set(z.hi.v,mpfr_cmp(l.v,r.v)>0?l.v:r.v,MPFR_RNDU);
    }
}
inline std::array<double,2> out(const Range& x) {
    return {mpfr_get_d(x.lo.v,MPFR_RNDD),mpfr_get_d(x.hi.v,MPFR_RNDU)};
}
} // namespace sproduct_detail

[[nodiscard]] inline SProductIntervals s_primitive_product_intervals(
    const SPrimitive& a,const SPrimitive& b) {
    if(!(a.exponent>0) || !(b.exponent>0) ||
       !std::isfinite(a.exponent) || !std::isfinite(b.exponent) ||
       !std::isfinite(a.coefficient) || !std::isfinite(b.coefficient))
        throw std::invalid_argument("invalid s primitive");
    for(int d=0;d<3;++d)
        if(!std::isfinite(a.center[d]) || !std::isfinite(b.center[d]))
            throw std::invalid_argument("invalid primitive center");
    using namespace sproduct_detail;
    Range alpha(a.exponent),beta(b.exponent),p,ab,mu;
    add(p,alpha,beta);
    mul(ab,alpha,beta);
    div_positive(mu,ab,p);
    Range distance2(0.0);
    SProductIntervals result{};
    for(int d=0;d<3;++d) {
        Range A(a.center[d]),B(b.center[d]),delta,delta2,next;
        sub(delta,A,B);
        square(delta2,delta);
        add(next,distance2,delta2);
        mpfr_set(distance2.lo.v,next.lo.v,MPFR_RNDD);
        mpfr_set(distance2.hi.v,next.hi.v,MPFR_RNDU);
        Range x,y,sum,center;
        mul(x,alpha,A);
        mul(y,beta,B);
        add(sum,x,y);
        div_positive(center,sum,p);
        result.center[d]=out(center);
    }
    Range md2;
    mul(md2,mu,distance2);
    Range decay;
    mpfr_neg(decay.lo.v,md2.hi.v,MPFR_RNDD);
    mpfr_neg(decay.hi.v,md2.lo.v,MPFR_RNDU);
    Range exponential;
    mpfr_exp(exponential.lo.v,decay.lo.v,MPFR_RNDD);
    mpfr_exp(exponential.hi.v,decay.hi.v,MPFR_RNDU);
    Range ca(std::abs(a.coefficient)),cb(std::abs(b.coefficient)),coeff,amp;
    mul(coeff,ca,cb);
    mul(amp,coeff,exponential);
    result.exponent=out(p);
    result.reduced_exponent=out(mu);
    result.amplitude=out(amp);
    // MPFR itself has a finite exponent range. A positive real product
    // must never acquire a zero upper bound solely from MPFR underflow.
    if (a.coefficient != 0.0 && b.coefficient != 0.0 &&
        result.amplitude[1] == 0.0)
        result.amplitude[1] = std::numeric_limits<double>::denorm_min();
    return result;
}
} // namespace chem::kumo::certified
