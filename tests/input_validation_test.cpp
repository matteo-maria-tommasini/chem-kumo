#include <chem/kumo.hpp>
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace chem::kumo;
template<class F> void expect_invalid(F&& f, const char* name) {
    try { f(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error(std::string("expected invalid_argument: ") + name);
}
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    try {
        auto mol = Molecule::from_fchk(argv[1]);
        auto g = Grid::around(mol.atoms(), 2.0, 0.5);
        expect_invalid([&]{ auto x = mol.densities_blocked(g,{Spin::Total},0); (void)x; }, "zero block");
        expect_invalid([&]{ auto x = mol.densities_blocked(g,{Spin::Total},65537); (void)x; }, "large block");
        auto bad = g; bad.counts[0] = 0;
        expect_invalid([&]{ auto x = mol.density(bad); (void)x; }, "zero grid");
        bad = g; bad.origin[0] = std::numeric_limits<double>::quiet_NaN();
        expect_invalid([&]{ auto x = mol.density(bad); (void)x; }, "NaN grid");
        bad = g; bad.steps.setZero();
        expect_invalid([&]{ auto x = mol.density(bad); (void)x; }, "singular grid");
        expect_invalid([&]{ auto x = Grid::around(mol.atoms(),-1.0,0.5); (void)x; }, "negative margin");
        expect_invalid([&]{ auto x = Grid::around(mol.atoms(),1.0,0.0); (void)x; }, "zero spacing");
        auto empty = mol.densities(g,{});
        if (!empty.empty()) throw std::runtime_error("empty request returned fields");
        auto dedup = mol.densities(g,{Spin::Total,Spin::Total});
        if (dedup.size()!=1) throw std::runtime_error("duplicate request not deduplicated");
        std::cout << "input validation OK\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
