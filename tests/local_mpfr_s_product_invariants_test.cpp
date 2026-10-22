#include "chem/kumo/local_mpfr_s_product.hpp"
#include <mpfr.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using chem::kumo::certified::SPrimitive;
using chem::kumo::certified::s_primitive_product_intervals;
struct M {
    mpfr_t v;
    explicit M(double x=0) { mpfr_init2(v,1024); mpfr_set_d(v,x,MPFR_RNDN); }
    ~M(){mpfr_clear(v);} M(const M&)=delete; M& operator=(const M&)=delete;
};
static void check_interval(const std::array<double,2>& bounds, mpfr_srcptr lo,
                           mpfr_srcptr hi, const char* label, int case_id) {
    if(std::isnan(bounds[0]) || std::isnan(bounds[1]) || bounds[0]>bounds[1])
        throw std::runtime_error(std::string(label)+" malformed case="+std::to_string(case_id));
    M lower(bounds[0]),upper(bounds[1]);
    if(mpfr_cmp(lower.v,lo)>0 || mpfr_cmp(upper.v,hi)<0)
        throw std::runtime_error(std::string(label)+" enclosure case="+std::to_string(case_id));
}
static void reference(const SPrimitive& a,const SPrimitive& b,
                      const chem::kumo::certified::SProductIntervals& out,int id,
                      std::size_t& checks) {
    M alpha(a.exponent),beta(b.exponent),pL,pU,muL,muU,prod,distL,distU;
    mpfr_add(pL.v,alpha.v,beta.v,MPFR_RNDD);
    mpfr_add(pU.v,alpha.v,beta.v,MPFR_RNDU);
    check_interval(out.exponent,pL.v,pU.v,"p",id); ++checks;
    mpfr_mul(prod.v,alpha.v,beta.v,MPFR_RNDD);
    mpfr_div(muL.v,prod.v,pU.v,MPFR_RNDD);
    mpfr_mul(prod.v,alpha.v,beta.v,MPFR_RNDU);
    mpfr_div(muU.v,prod.v,pL.v,MPFR_RNDU);
    check_interval(out.reduced_exponent,muL.v,muU.v,"mu",id); ++checks;
    mpfr_set_zero(distL.v,0);mpfr_set_zero(distU.v,0);
    for(int d=0;d<3;++d) {
        M A(a.center[d]),B(b.center[d]),x,y,z,cl,cu,deltaL,deltaU,sqL,sqU;
        mpfr_mul(x.v,alpha.v,A.v,MPFR_RNDD);
        mpfr_mul(y.v,beta.v,B.v,MPFR_RNDD);
        mpfr_add(z.v,x.v,y.v,MPFR_RNDD);
        mpfr_div(cl.v,z.v,pU.v,MPFR_RNDD);
        mpfr_mul(x.v,alpha.v,A.v,MPFR_RNDU);
        mpfr_mul(y.v,beta.v,B.v,MPFR_RNDU);
        mpfr_add(z.v,x.v,y.v,MPFR_RNDU);
        mpfr_div(cu.v,z.v,pL.v,MPFR_RNDU);
        // For signed numerators the quotient extrema require interval division.
        // Recompute using all four endpoint quotients, with directed rounding.
        M nL,nU,t;
        mpfr_mul(x.v,alpha.v,A.v,MPFR_RNDD);mpfr_mul(y.v,beta.v,B.v,MPFR_RNDD);
        mpfr_add(nL.v,x.v,y.v,MPFR_RNDD);
        mpfr_mul(x.v,alpha.v,A.v,MPFR_RNDU);mpfr_mul(y.v,beta.v,B.v,MPFR_RNDU);
        mpfr_add(nU.v,x.v,y.v,MPFR_RNDU);
        mpfr_div(cl.v,nL.v,pL.v,MPFR_RNDD);
        mpfr_div(t.v,nL.v,pU.v,MPFR_RNDD);if(mpfr_cmp(t.v,cl.v)<0)mpfr_set(cl.v,t.v,MPFR_RNDD);
        mpfr_div(cu.v,nU.v,pL.v,MPFR_RNDU);
        mpfr_div(t.v,nU.v,pU.v,MPFR_RNDU);if(mpfr_cmp(t.v,cu.v)>0)mpfr_set(cu.v,t.v,MPFR_RNDU);
        check_interval(out.center[d],cl.v,cu.v,"center",id); ++checks;
        mpfr_sub(deltaL.v,A.v,B.v,MPFR_RNDD);
        mpfr_sub(deltaU.v,A.v,B.v,MPFR_RNDU);
        mpfr_sqr(sqL.v,deltaL.v,MPFR_RNDD);
        mpfr_sqr(sqU.v,deltaU.v,MPFR_RNDU);
        mpfr_add(distL.v,distL.v,sqL.v,MPFR_RNDD);
        mpfr_add(distU.v,distU.v,sqU.v,MPFR_RNDU);
    }
    M coeffA(std::abs(a.coefficient)),coeffB(std::abs(b.coefficient));
    M coeffL,coeffU,decL,decU,eL,eU,ampL,ampU;
    mpfr_mul(coeffL.v,coeffA.v,coeffB.v,MPFR_RNDD);
    mpfr_mul(coeffU.v,coeffA.v,coeffB.v,MPFR_RNDU);
    mpfr_mul(decL.v,muU.v,distU.v,MPFR_RNDU);mpfr_neg(decL.v,decL.v,MPFR_RNDD);
    mpfr_mul(decU.v,muL.v,distL.v,MPFR_RNDD);mpfr_neg(decU.v,decU.v,MPFR_RNDU);
    mpfr_exp(eL.v,decL.v,MPFR_RNDD);mpfr_exp(eU.v,decU.v,MPFR_RNDU);
    mpfr_mul(ampL.v,coeffL.v,eL.v,MPFR_RNDD);
    mpfr_mul(ampU.v,coeffU.v,eU.v,MPFR_RNDU);
    check_interval(out.amplitude,ampL.v,ampU.v,"amplitude",id); ++checks;
}
int main(){
 try {
    std::uint64_t state=0x9e3779b97f4a7c15ULL;
    auto rnd=[&](){state^=state>>12;state^=state<<25;state^=state>>27;
        return (state*2685821657736338717ULL)>>11;};
    auto unit=[&](){return static_cast<double>(rnd())*0x1.0p-53;};
    std::size_t checks=0;
    for(int i=0;i<1200;++i){
        auto expo=[&](){return std::pow(10.,-12.+24.*unit());};
        auto coord=[&](){return (2.*unit()-1.)*std::pow(10.,-4.+8.*unit());};
        auto coeff=[&](){return (2.*unit()-1.)*std::pow(10.,-5.+10.*unit());};
        SPrimitive a{expo(),coeff(),{coord(),coord(),coord()}};
        SPrimitive b{expo(),coeff(),{coord(),coord(),coord()}};
        if(i%11==0)b.center=a.center;
        if(i%17==0)b.coefficient=0.;
        if(i%29==0){a.exponent=std::numeric_limits<double>::denorm_min();}
        const auto out=s_primitive_product_intervals(a,b);
        reference(a,b,out,i,checks);
    }
    std::cout<<"cases=1200 checks="<<checks<<" reference_precision=1024 status=OK\n";
 }catch(const std::exception& e){std::cerr<<"status=FAIL reason="<<e.what()<<'\n';return 1;}
}
