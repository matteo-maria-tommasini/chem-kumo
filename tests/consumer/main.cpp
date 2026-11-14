#include <chem/kumo.hpp>
#include <type_traits>

int main() {
    using chem::kumo::DensityComponent;
    static_assert(std::is_enum_v<DensityComponent>);
    chem::kumo::Grid grid;
    return grid.counts.x() == 1 ? 0 : 1;
}
