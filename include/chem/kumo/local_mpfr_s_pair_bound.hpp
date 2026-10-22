#pragma once
// Certified overlap bound for degree-zero Gaussian pair intermediates.
// Certification is relative to the exact real values of the stored binary64
// fields, NOT to the original unrounded Gaussian primitives.
#include "chem/kumo/local_gaussian_pair.hpp"
#include <mpfr.h>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace chem::kumo::certified {
namespace detail {
struct Real {
    mpfr_t x;
    Real() { mpfr_init2(x, 256); }
    explicit Real(double v) : Real() { mpfr_set_d(x, v, MPFR_RNDN); }
    ~Real() { mpfr_clear(x); }
    Real(const Real&) = delete;
    Real& operator=(const Real&) = delete;
};
struct Range {
    Real lo, hi;
    explicit Range(double v) { mpfr_set_d(lo.x,v,MPFR_RNDN); mpfr_set_d(hi.x,v,MPFR_RNDN); }
    Range() = default;
    Range(const Range&) = delete;
    Range& operator=(const Range&) = delete;
};
inline void add(Range& z, const Range& a, const Range& b) {
    mpfr_add(z.lo.x,a.lo.x,b.lo.x,MPFR_RNDD);
    mpfr_add(z.hi.x,a.hi.x,b.hi.x,MPFR_RNDU);
}
inline void sub(Range& z, const Range& a, const Range& b) {
    mpfr_sub(z.lo.x,a.lo.x,b.hi.x,MPFR_RNDD);
    mpfr_sub(z.hi.x,a.hi.x,b.lo.x,MPFR_RNDU);
}
inline void mul_nonnegative(Range& z, const Range& a, const Range& b) {
    mpfr_mul(z.lo.x,a.lo.x,b.lo.x,MPFR_RNDD);
    mpfr_mul(z.hi.x,a.hi.x,b.hi.x,MPFR_RNDU);
}
inline void div_positive(Range& z, const Range& a, const Range& b) {
    if (mpfr_sgn(b.lo.x) <= 0) throw std::domain_error("nonpositive divisor interval");
    mpfr_div(z.lo.x,a.lo.x,b.hi.x,MPFR_RNDD);
    mpfr_div(z.hi.x,a.hi.x,b.lo.x,MPFR_RNDU);
}
} // namespace detail

// For f(r)=c_a exp(-a|r-A|^2) and g(r)=c_b exp(-b|r-B|^2),
// returns a rounded-up upper bound on integral |f(r)g(r)| dr.
// All stored input fields are interpreted as exact binary64 numbers.
[[nodiscard]] inline double s_pair_overlap_upper(
    const GaussianPairIntermediate& a, const GaussianPairIntermediate& b) {
    for (const auto* p : {&a, &b}) {
        if (!(p->exponent > 0) || !std::isfinite(p->exponent) ||
            !std::isfinite(p->prefactor))
            throw std::invalid_argument("invalid Gaussian pair parameters");
        for (int d=0;d<3;++d)
            if (p->degree[d] != 0 || !std::isfinite(p->center[d]) ||
                !std::isfinite(p->polynomial[d][0]))
                throw std::invalid_argument("requires finite degree-zero pair");
    }
    using namespace detail;
    Range ea(a.exponent), eb(b.exponent), p, ab, eta;
    add(p,ea,eb);
    mul_nonnegative(ab,ea,eb);
    div_positive(eta,ab,p);
    Range distance2(0.0);
    for (int d=0;d<3;++d) {
        Range x(a.center[d]), y(b.center[d]), delta, square;
        sub(delta,x,y);
        // Square an interval without introducing a spurious negative lower bound.
        if (mpfr_sgn(delta.lo.x) >= 0) {
            mpfr_sqr(square.lo.x,delta.lo.x,MPFR_RNDD);
            mpfr_sqr(square.hi.x,delta.hi.x,MPFR_RNDU);
        } else if (mpfr_sgn(delta.hi.x) <= 0) {
            mpfr_sqr(square.lo.x,delta.hi.x,MPFR_RNDD);
            mpfr_sqr(square.hi.x,delta.lo.x,MPFR_RNDU);
        } else {
            mpfr_set_zero(square.lo.x,0);
            Real l,r;
            mpfr_sqr(l.x,delta.lo.x,MPFR_RNDU);
            mpfr_sqr(r.x,delta.hi.x,MPFR_RNDU);
            mpfr_set(square.hi.x,mpfr_cmp(l.x,r.x)>0?l.x:r.x,MPFR_RNDU);
        }
        Range next;
        add(next,distance2,square);
        mpfr_set(distance2.lo.x,next.lo.x,MPFR_RNDD);
        mpfr_set(distance2.hi.x,next.hi.x,MPFR_RNDU);
    }
    Range decay_exponent;
    mul_nonnegative(decay_exponent,eta,distance2);
    Real decay_upper;
    mpfr_neg(decay_upper.x,decay_exponent.lo.x,MPFR_RNDU);
    mpfr_exp(decay_upper.x,decay_upper.x,MPFR_RNDU);
    Real pi, ratio, root, volume;
    mpfr_const_pi(pi.x,MPFR_RNDU);
    mpfr_div(ratio.x,pi.x,p.lo.x,MPFR_RNDU);
    mpfr_sqrt(root.x,ratio.x,MPFR_RNDU);
    mpfr_mul(volume.x,ratio.x,root.x,MPFR_RNDU);
    Real result;
    mpfr_set_d(result.x,std::abs(a.prefactor),MPFR_RNDU);
    mpfr_mul(result.x,result.x,decay_upper.x,MPFR_RNDU);
    mpfr_mul(result.x,result.x,volume.x,MPFR_RNDU);
    mpfr_mul_d(result.x,result.x,std::abs(b.prefactor),MPFR_RNDU);
    for (int d=0;d<3;++d) {
        mpfr_mul_d(result.x,result.x,std::abs(a.polynomial[d][0]),MPFR_RNDU);
        mpfr_mul_d(result.x,result.x,std::abs(b.polynomial[d][0]),MPFR_RNDU);
    }
    const double upper=mpfr_get_d(result.x,MPFR_RNDU);
    if (std::isnan(upper) || upper<0) throw std::runtime_error("invalid certified upper bound");
    return upper;
}
} // namespace chem::kumo::certified
