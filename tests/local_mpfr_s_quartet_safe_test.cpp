#include "chem/kumo/local_mpfr_s_quartet_safe.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace chem::kumo;
using namespace chem::kumo::certified;
static void check(bool condition, const char* what) {
    if (!condition) throw std::runtime_error(what);
}
int main() {
    try {
        const SPrimitive a{1.0,1.0,{0,0,0}};
        const SPrimitive b{2.0,1.0,{1,0,0}};
        const auto good=safe_s_primitive_quartet_bound(a,b,a,b);
        check(good.status==SQuartetBoundStatus::FiniteEngineering,"finite status");
        check(std::isfinite(good.bound.value) && good.bound.value>0,"finite upper");
        check(!local_bound_allows_skip(good.bound,1e100,BoundPolicy::Rigorous),"rigorous gate");
        check(local_bound_allows_skip(good.bound,1e100,BoundPolicy::Experimental),"experimental gate");
        auto bad=a; bad.exponent=std::numeric_limits<double>::quiet_NaN();
        const auto invalid=safe_s_primitive_quartet_bound(bad,b,a,b);
        check(invalid.status==SQuartetBoundStatus::InvalidInput,"invalid status");
        check(!local_bound_allows_skip(invalid.bound,1e100,BoundPolicy::Experimental),"invalid gate");
        auto zero=a; zero.coefficient=0;
        const auto exact=safe_s_primitive_quartet_bound(zero,b,a,b);
        check(exact.status==SQuartetBoundStatus::ExactZeroCoefficient,"exact zero status");
        check(exact.bound.value==0,"exact zero value");
        check(!local_bound_allows_skip(exact.bound,1e100,BoundPolicy::Experimental),"zero gate");
        auto far=a; far.center={1e100,0,0};
        const auto extreme=safe_s_primitive_quartet_bound(a,b,far,b);
        check(extreme.status==SQuartetBoundStatus::FiniteEngineering ||
              extreme.status==SQuartetBoundStatus::ConservativeFallback,"extreme status");
        check(!local_bound_allows_skip(extreme.bound,1e100,BoundPolicy::Rigorous),"extreme rigorous gate");
        std::cout << "checks=10 fail_closed=OK rigorous=keep experimental=OK status=OK\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
