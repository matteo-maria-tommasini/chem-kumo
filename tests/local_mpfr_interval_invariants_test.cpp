#include "chem/kumo/local_mpfr_s_product.hpp"
#include <mpfr.h>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstdint>
using namespace chem::kumo::certified::sproduct_detail;

struct Ref {
    mpfr_t v;
    explicit Ref(double x=0.0) { mpfr_init2(v,1024); mpfr_set_d(v,x,MPFR_RNDN); }
    ~Ref() { mpfr_clear(v); }
    Ref(const Ref&)=delete;
    Ref& operator=(const Ref&)=delete;
};
static void require(bool ok, const char* label, std::size_t i) {
    if(!ok) throw std::runtime_error(std::string(label)+" case="+std::to_string(i));
}
static bool encloses(const Range& x, mpfr_srcptr exact) {
    return !mpfr_nan_p(x.lo.v) && !mpfr_nan_p(x.hi.v) &&
           mpfr_cmp(x.lo.v,x.hi.v)<=0 &&
           mpfr_cmp(x.lo.v,exact)<=0 && mpfr_cmp(x.hi.v,exact)>=0;
}
static void verify(const char* name, const Range& x, mpfr_srcptr exact, std::size_t i,
                   std::size_t& checks, std::size_t& nonfinite) {
    require(encloses(x,exact),name,i);
    if(!mpfr_number_p(x.lo.v) || !mpfr_number_p(x.hi.v)) ++nonfinite;
    ++checks;
}
int main() {
 try {
    const std::array<double,19> v = {
        0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 3.0, -3.0,
        std::numeric_limits<double>::denorm_min(),
        -std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::min(),
        -std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(),
        1e-100, -1e-100, 1e100, -1e100, 1e-300
    };
    std::size_t checks=0, nonfinite=0, cases=0;
    Ref exact;
    for(double a:v) for(double b:v) {
        const std::size_t i=cases++;
        Range x(a), y(b), z;
        Ref aa(a),bb(b);
        add(z,x,y); mpfr_add(exact.v,aa.v,bb.v,MPFR_RNDN);
        verify("add",z,exact.v,i,checks,nonfinite);
        sub(z,x,y); mpfr_sub(exact.v,aa.v,bb.v,MPFR_RNDN);
        verify("sub",z,exact.v,i,checks,nonfinite);
        mul(z,x,y); mpfr_mul(exact.v,aa.v,bb.v,MPFR_RNDN);
        verify("mul",z,exact.v,i,checks,nonfinite);
        square(z,x); mpfr_sqr(exact.v,aa.v,MPFR_RNDN);
        verify("square",z,exact.v,i,checks,nonfinite);
        if(b>0) {
            div_positive(z,x,y); mpfr_div(exact.v,aa.v,bb.v,MPFR_RNDN);
            verify("div_positive",z,exact.v,i,checks,nonfinite);
        }
    }
    // Explicit crossing-zero square interval, including a very asymmetric one.
    for(const auto& pair: std::array<std::array<double,2>,3>{{{-3.0,2.0},{-1e100,1e-100},{-1e-300,1e100}}}) {
        Range x,z; mpfr_set_d(x.lo.v,pair[0],MPFR_RNDN);
        mpfr_set_d(x.hi.v,pair[1],MPFR_RNDN);
        square(z,x);
        require(mpfr_zero_p(z.lo.v),"cross-zero square lower",cases);
        Ref l(pair[0]),r(pair[1]);
        mpfr_sqr(l.v,l.v,MPFR_RNDN); mpfr_sqr(r.v,r.v,MPFR_RNDN);
        require(mpfr_cmp(z.hi.v,l.v)>=0 && mpfr_cmp(z.hi.v,r.v)>=0,
                "cross-zero square upper",cases);
        ++checks; ++cases;
    }
    std::cout<<"cases="<<cases<<" checks="<<checks
             <<" nonfinite_intervals="<<nonfinite
             <<" reference_precision=1024 status=OK\n";
    return 0;
 } catch(const std::exception& e) {
    std::cerr<<"status=FAIL reason="<<e.what()<<'\n'; return 1;
 }
}
