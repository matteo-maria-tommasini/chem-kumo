#include "backend_selection.hpp"
#include "chem/kumo/density.hpp"
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <stdexcept>
using namespace chem::kumo::detail;
int main(int argc, char** argv) {
    if (argc != 3) return 10;
    for (int i = 1; i < argc; ++i) {
        const auto molecule = chem::kumo::Molecule::from_fchk(argv[i]);
        const auto& info = molecule.electronic_structure_info();
        if (info.occupied_alpha != static_cast<std::size_t>(molecule.occupied_count(chem::kumo::Spin::Alpha)) ||
            info.occupied_beta != static_cast<std::size_t>(molecule.occupied_count(chem::kumo::Spin::Beta))) return 11;
        if (!std::isfinite(info.total_relative_error) || !std::isfinite(info.spin_relative_error)) return 12;
        if (info.ao_mo_eligible() != (info.orbital_status == chem::kumo::OrbitalRepresentationStatus::Valid)) return 13;
        std::cout << argv[i] << " status=" << static_cast<int>(info.orbital_status)
                  << " total_abs=" << info.total_absolute_error
                  << " total_rel=" << info.total_relative_error
                  << " spin_abs=" << info.spin_absolute_error
                  << " spin_rel=" << info.spin_relative_error << '\n';
        const auto second = chem::kumo::Molecule::from_fchk(argv[i]);
        if (second.electronic_structure_info().orbital_status != info.orbital_status) return 14;
    }
    Eigen::MatrixXd a = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd b = Eigen::MatrixXd::Identity(3, 3);
    b.col(0) = Eigen::Vector3d(0, 1, 0);
    Eigen::MatrixXd pa = a.leftCols(1) * a.leftCols(1).transpose();
    Eigen::MatrixXd pb = b.leftCols(1) * b.leftCols(1).transpose();
    auto yes = validate_orbital_backend(pa + pb, pa - pb, a, b, 1, 1, true);
    if (yes.backend != ElectronicBackend::AoMo) return 1;
    auto no = validate_orbital_backend(pa + pb, pa - pb + Eigen::MatrixXd::Identity(3, 3)*1e-4, a, b, 1, 1, true);
    if (no.backend != ElectronicBackend::AoP) return 2;
    auto restricted = validate_orbital_backend(2*pa, Eigen::MatrixXd::Zero(3, 3), a, Eigen::MatrixXd{}, 1, 1, false);
    if (restricted.backend != ElectronicBackend::AoMo) return 3;
    if (validate_orbital_backend(pa+pb, pa-pb, a, b, 4, 1, true).backend != ElectronicBackend::AoP) return 4;
    if (validate_orbital_backend(pa+pb, pa-pb, a, b, 1, 1, false).backend != ElectronicBackend::AoP) return 5;
    try { (void)validate_orbital_backend(pa+pb, pa-pb, a, b, 1, 1, true, -1); }
    catch (const std::invalid_argument&) { return 0; }
    return 6;
}
