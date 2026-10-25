#include "chem/kumo/field_integrals.hpp"
#include "chem/kumo/field_sampler.hpp"
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

int main() {
    chem::kumo::Grid grid;
    grid.counts = {2, 2, 2};
    grid.steps << 1, .2, 0, 0, 1, 0, 0, 0, 1;
    chem::kumo::Field field(grid, {}, {1,2,3,4,5,6,7,8});
    assert(std::abs(chem::kumo::grid_cell_volume(field)-1) < 1e-14);
    assert(std::abs(chem::kumo::integrate_rectangular(field)-36) < 1e-14);
    assert(std::abs(chem::kumo::integrate_trapezoidal(field)-4.5) < 1e-14);
    grid.counts = {1, 2, 2};
    chem::kumo::Field single(grid, {}, {1,2,3,4});
    assert(std::abs(chem::kumo::integrate_rectangular(single)-10) < 1e-14);
    assert(std::abs(chem::kumo::integrate_trapezoidal(single)-2.5) < 1e-14);
    bool rejected = false;
    try { (void)chem::kumo::FieldSampler(single); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    single.mutable_values()[0] = std::numeric_limits<double>::quiet_NaN();
    rejected = false;
    try { (void)chem::kumo::integrate_trapezoidal(single); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}
