#include <chem/kumo.hpp>
#include <iostream>
int main() {
    chem::kumo::Grid grid;
    if (grid.counts.x() != 1) return 1;
    std::cout << "installed KUMO API OK\n";
    return 0;
}
