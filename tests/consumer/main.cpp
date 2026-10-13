#include <chem/kumo.hpp>
#include <type_traits>
int main() {
    using chem::kumo::Spin;
    static_assert(std::is_enum_v<Spin>);
    chem::kumo::Grid grid;
    return grid.counts.x() == 1 ? 0 : 1;
}
