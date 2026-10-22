#include "chem/kumo/local_sparse_budget.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using chem::kumo::detail::LocalBudgetCharge;
using chem::kumo::detail::LocalSparseBudget;

int main() {
    try {
        int checks = 0;
        auto check = [&](bool ok) {
            ++checks;
            if (!ok) throw std::runtime_error("extreme budget regression failed at check " + std::to_string(checks));
        };
        const double tiny = std::numeric_limits<double>::denorm_min();
        const double next = std::nextafter(1.0, std::numeric_limits<double>::infinity());
        LocalSparseBudget zero(4, {0.0, 0.0});
        check(zero.try_charge({}));
        check(zero.try_charge({{0,0,0,0.0},{1,2,1,0.0}}));
        check(!zero.try_charge({{0,0,0,tiny}}));
        check(zero.spent()[0](0,0) == 0.0);
        check(!zero.try_charge({{0,0,0,tiny},{1,2,1,tiny}}));
        check(zero.spent()[1](1,2) == 0.0);
        LocalSparseBudget one(4, {1.0, tiny});
        check(one.try_charge({{0,0,0,1.0}}));
        check(!one.try_charge({{0,0,0,next-1.0}}));
        check(one.spent()[0](0,0) == 1.0);
        check(one.try_charge({{0,0,1,tiny}}));
        check(!one.try_charge({{0,0,1,tiny}}));
        check(one.spent()[1](0,0) == tiny);
        LocalSparseBudget duplicates(4, {1.0});
        check(!duplicates.try_charge({{0,0,0,0.6},{0,0,0,0.6}}));
        check(duplicates.spent()[0](0,0) == 0.0);
        check(duplicates.try_charge({{0,0,0,0.5},{0,0,0,0.5}}));
        check(!duplicates.try_charge({{0,0,0,0.25},{0,1,0,tiny}}));
        check(duplicates.spent()[0](0,1) == 0.0);
        LocalSparseBudget near(4, {next});
        check(near.try_charge({{0,0,0,1.0}}));
        check(near.try_charge({{0,0,0,next-1.0}}));
        check(near.spent()[0](0,0) == next);
        check(!near.try_charge({{0,0,0,next-1.0}}));
        auto invalid = [&](double x) {
            try { one.try_charge({{0,0,0,x}}); return false; }
            catch (const std::invalid_argument&) { return true; }
        };
        check(invalid(std::numeric_limits<double>::quiet_NaN()));
        check(invalid(std::numeric_limits<double>::infinity()));
        check(invalid(-tiny));
        std::cout << "checks=" << checks << " zero_budget=OK subnormal=OK duplicate=OK atomicity=OK invalid=OK status=OK\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
