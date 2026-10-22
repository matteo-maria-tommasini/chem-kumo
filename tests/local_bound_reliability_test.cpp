#include "chem/kumo/local_bound_reliability.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace chem::kumo;

static void check(bool value) {
    if (!value) throw std::runtime_error("bound reliability regression failed");
}
int main() {
    try {
        const double nan=std::numeric_limits<double>::quiet_NaN();
        const double inf=std::numeric_limits<double>::infinity();
        const double sub=std::numeric_limits<double>::denorm_min();
        const ReliableBound engineering{0.25,BoundReliability::Engineering};
        const ReliableBound certified{0.25,BoundReliability::Certified};
        check(!local_bound_allows_skip(engineering,1.0));
        check(local_bound_allows_skip(engineering,1.0,BoundPolicy::Experimental));
        check(local_bound_allows_skip(certified,1.0));
        check(!local_bound_allows_skip(certified,0.2));
        check(!local_bound_allows_skip(certified,0.0));
        check(!local_bound_allows_skip(certified,-1.0));
        check(!local_bound_allows_skip(certified,nan));
        check(!local_bound_allows_skip(certified,inf));
        check(!local_bound_allows_skip({nan,BoundReliability::Certified},1.0));
        check(!local_bound_allows_skip({inf,BoundReliability::Certified},1.0));
        check(!local_bound_allows_skip({-1.0,BoundReliability::Certified},1.0));
        check(!local_bound_allows_skip({0.0,BoundReliability::Certified},1.0));
        check(!local_bound_allows_skip({0.0,BoundReliability::Engineering},1.0,BoundPolicy::Experimental));
        check(!local_bound_allows_skip({sub,BoundReliability::Unusable},1.0));
        check(!local_bound_allows_skip({0.1,BoundReliability::Unusable},1.0,BoundPolicy::Experimental));
        check(local_bound_allows_skip({sub,BoundReliability::Certified},sub));
        check(!local_bound_allows_skip({sub,BoundReliability::Engineering},sub));
        check(local_bound_allows_skip({sub,BoundReliability::Engineering},sub,BoundPolicy::Experimental));
        std::cout << "checks=18 rigorous=OK experimental=OK zero=keep invalid=keep status=OK\n";
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
