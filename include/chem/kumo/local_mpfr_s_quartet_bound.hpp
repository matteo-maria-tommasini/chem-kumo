#pragma once
// Experimental MPFR outward-rounded overlap bound for four s primitives.
// Input binary64 values are interpreted as exact real numbers.
#include "chem/kumo/local_mpfr_s_product.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace chem::kumo::certified {
namespace squartet_detail {
using namespace sproduct_detail;

inline void copy(Range& dst, const Range& src) {
    mpfr_set(dst.lo.v,src.lo.v,MPFR_RNDD);
    mpfr_set(dst.hi.v,src.hi.v,MPFR_RNDU);
}
inline void positive_product(Range& z, const Range& x, const Range& y) {
    mpfr_mul(z.lo.v,x.lo.v,y.lo.v,MPFR_RNDD);
    mpfr_mul(z.hi.v,x.hi.v,y.hi.v,MPFR_RNDU);
}
inline void product_data(const SPrimitive& a, const SPrimitive& b,
                         Range& p, Range& mu, Range (&center)[3], Range& d2) {
    Range alpha(a.exponent),beta(b.exponent),ab;
    add(p,alpha,beta);
    positive_product(ab,alpha,beta);
    div_positive(mu,ab,p);
    for(int d=0;d<3;++d) {
        Range A(a.center[d]),B(b.center[d]),diff,dsq,tmp;
        sub(diff,A,B);
        square(dsq,diff);
        add(tmp,d2,dsq);
        copy(d2,tmp);
        Range aa,bb,sum;
        mul(aa,alpha,A);
        mul(bb,beta,B);
        add(sum,aa,bb);
        div_positive(center[d],sum,p);
    }
}
} // namespace squartet_detail

[[nodiscard]] inline double s_primitive_quartet_overlap_upper(
    const SPrimitive& a,const SPrimitive& b,
    const SPrimitive& c,const SPrimitive& d) {
    const SPrimitive* primitives[]={&a,&b,&c,&d};
    for(const auto* v:primitives) {
        if(!(v->exponent>0.0) || !std::isfinite(v->exponent) ||
           !std::isfinite(v->coefficient))
            throw std::invalid_argument("invalid s primitive exponent or coefficient");
        for(double x:v->center)
            if(!std::isfinite(x)) throw std::invalid_argument("invalid s primitive center");
    }
    if(a.coefficient==0 || b.coefficient==0 || c.coefficient==0 || d.coefficient==0)
        return 0.0;
    using namespace squartet_detail;
    Range p,mu,pa[3],dab2(0.0),q,nu,pb[3],dcd2(0.0);
    product_data(a,b,p,mu,pa,dab2);
    product_data(c,d,q,nu,pb,dcd2);
    Range total, pq, eta, pq2(0.0);
    add(total,p,q);
    positive_product(pq,p,q);
    div_positive(eta,pq,total);
    for(int axis=0;axis<3;++axis) {
        Range diff,dsq,tmp;
        sub(diff,pa[axis],pb[axis]);
        square(dsq,diff);
        add(tmp,pq2,dsq);
        copy(pq2,tmp);
    }
    Range t1,t2,t3,sum12,decay;
    positive_product(t1,mu,dab2);
    positive_product(t2,nu,dcd2);
    positive_product(t3,eta,pq2);
    add(sum12,t1,t2);
    add(decay,sum12,t3);
    // Work in log space: exp(-decay) may underflow MPFR even when
    // the final upward-rounded binary64 bound is simply denorm_min().
    // Every contribution to log_upper is rounded towards +infinity.
    Num exponent,pi,ratio,log_ratio,log_volume,log_upper,term,answer;
    mpfr_neg(exponent.v,decay.lo.v,MPFR_RNDU);
    mpfr_const_pi(pi.v,MPFR_RNDU);
    mpfr_div(ratio.v,pi.v,total.lo.v,MPFR_RNDU);
    mpfr_log(log_ratio.v,ratio.v,MPFR_RNDU);
    mpfr_mul_d(log_volume.v,log_ratio.v,1.5,MPFR_RNDU);
    mpfr_add(log_upper.v,exponent.v,log_volume.v,MPFR_RNDU);
    for(const auto* v:primitives) {
        mpfr_set_d(term.v,std::abs(v->coefficient),MPFR_RNDN);
        mpfr_log(term.v,term.v,MPFR_RNDU);
        mpfr_add(log_upper.v,log_upper.v,term.v,MPFR_RNDU);
    }
    if(mpfr_nan_p(log_upper.v))
        return std::numeric_limits<double>::infinity();
    // Directed bounds for the binary64 representability thresholds.
    // If log_upper <= floor(log(denorm_min)), denorm_min is a safe
    // positive upper bound.  The lower-rounded logarithm is essential:
    // comparing with an upper-rounded threshold would be unsafe.
    Num log_min_subnormal, log_max_finite;
    mpfr_set_d(log_min_subnormal.v,
               std::numeric_limits<double>::denorm_min(),MPFR_RNDN);
    mpfr_log(log_min_subnormal.v,log_min_subnormal.v,MPFR_RNDD);
    if(mpfr_cmp(log_upper.v,log_min_subnormal.v)<=0)
        return std::numeric_limits<double>::denorm_min();

    // If log_upper > ceil(log(DBL_MAX)), the computed upper bound
    // cannot be represented by a finite binary64 value.  Returning
    // infinity is conservative; it says nothing about the exact integral.
    mpfr_set_d(log_max_finite.v,
               std::numeric_limits<double>::max(),MPFR_RNDN);
    mpfr_log(log_max_finite.v,log_max_finite.v,MPFR_RNDU);
    if(mpfr_cmp(log_upper.v,log_max_finite.v)>0)
        return std::numeric_limits<double>::infinity();
    if(!mpfr_number_p(log_upper.v))
        return std::numeric_limits<double>::infinity();
    mpfr_clear_flags();
    mpfr_exp(answer.v,log_upper.v,MPFR_RNDU);
    if(mpfr_nan_p(answer.v) || mpfr_overflow_p() || mpfr_underflow_p())
        return std::numeric_limits<double>::infinity();
    const double upper=mpfr_get_d(answer.v,MPFR_RNDU);
    if(upper==0.0) return std::numeric_limits<double>::denorm_min();
    return upper;
}
} // namespace chem::kumo::certified
