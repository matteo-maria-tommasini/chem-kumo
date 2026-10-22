#include "chem/kumo/local_mpfr_interval.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace chem::kumo::certified;
int main() {
    int checks=0;
    auto check=[&](bool b) { ++checks; if(!b) throw std::runtime_error("MPFR interval check failed"); };
    const auto add=evaluate(Operation::Add,0.1,0.2);
    check(add.lower <= (0.1+0.2) && add.upper >= (0.1+0.2));
    check(add.lower < add.upper);
    const auto mul=evaluate(Operation::Multiply,0.1,0.2);
    check(mul.lower <= mul.upper && mul.lower > 0.0);
    const auto div=evaluate(Operation::Divide,1.0,3.0);
    check(div.lower < div.upper && div.lower <= 1.0/3.0 && div.upper >= 1.0/3.0);
    const auto sq=evaluate(Operation::SquareRoot,2.0);
    check(sq.lower < sq.upper && sq.lower <= std::sqrt(2.0) && sq.upper >= std::sqrt(2.0));
    const auto ex=evaluate(Operation::Exponential,-1000.0);
    check(ex.lower == 0.0 && ex.upper > 0.0);
    const auto ex2=evaluate(Operation::Exponential,1000.0);
    check(ex2.lower == std::numeric_limits<double>::max() &&
          std::isinf(ex2.upper) && ex2.upper > 0.0);
    const auto exact=evaluate(Operation::Multiply,2.0,4.0);
    check(exact.lower == 8.0 && exact.upper == 8.0);
    bool rejected=false;
    try { (void)evaluate(Operation::Divide,1.0,0.0); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected);
    rejected=false;
    try { (void)evaluate(Operation::SquareRoot,-1.0); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected);
    std::cout << "checks=" << checks << " mpfr_directed=OK subnormal=OK overflow=OK status=OK\n";
}
