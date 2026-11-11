// AO-MO eligibility: the occupied MOs must reproduce the FCHK density
// matrices. Checked on real FCHK files and on small synthetic matrices.
#include "backend_selection.hpp"
#include "chem/kumo/density.hpp"
#include "test_support.hpp"
#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <stdexcept>
using namespace chem::kumo::detail;
using kumo_test::check;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " FILE1.fchk FILE2.fchk\n";
        return 2;
    }
    for (int i = 1; i < argc; ++i) {
        const std::string file = argv[i];
        const auto molecule = chem::kumo::Molecule::from_fchk(file);
        const auto& info = molecule.electronic_structure_info();
        check(info.occupied_alpha ==
                      static_cast<std::size_t>(molecule.occupied_count(chem::kumo::Spin::Alpha)) &&
                  info.occupied_beta ==
                      static_cast<std::size_t>(molecule.occupied_count(chem::kumo::Spin::Beta)),
              file + ": occupied counts in ElectronicStructureInfo match the molecule");
        check(std::isfinite(info.total_relative_error) && std::isfinite(info.spin_relative_error),
              file + ": representation errors are finite");
        check(info.ao_mo_eligible() ==
                  (info.orbital_status == chem::kumo::OrbitalRepresentationStatus::Valid),
              file + ": ao_mo_eligible() agrees with orbital_status");
        std::cout << file << " status=" << static_cast<int>(info.orbital_status)
                  << " total_abs=" << info.total_absolute_error
                  << " total_rel=" << info.total_relative_error
                  << " spin_abs=" << info.spin_absolute_error
                  << " spin_rel=" << info.spin_relative_error << '\n';
        const auto second = chem::kumo::Molecule::from_fchk(file);
        check(second.electronic_structure_info().orbital_status == info.orbital_status,
              file + ": eligibility is deterministic across loads");
    }

    // Synthetic unrestricted case: one alpha and one beta electron.
    Eigen::MatrixXd a = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd b = Eigen::MatrixXd::Identity(3, 3);
    b.col(0) = Eigen::Vector3d(0, 1, 0);
    Eigen::MatrixXd pa = a.leftCols(1) * a.leftCols(1).transpose();
    Eigen::MatrixXd pb = b.leftCols(1) * b.leftCols(1).transpose();
    check(validate_orbital_backend(pa + pb, pa - pb, a, b, 1, 1, true).backend ==
              ElectronicBackend::AoMo,
          "consistent unrestricted MOs are AO-MO eligible");
    check(validate_orbital_backend(pa + pb, pa - pb + Eigen::MatrixXd::Identity(3, 3) * 1e-4, a, b,
                                   1, 1, true)
                  .backend == ElectronicBackend::AoP,
          "a perturbed spin density falls back to AO-P");
    check(validate_orbital_backend(2 * pa, Eigen::MatrixXd::Zero(3, 3), a, Eigen::MatrixXd{}, 1, 1,
                                   false)
                  .backend == ElectronicBackend::AoMo,
          "consistent restricted MOs are AO-MO eligible");
    check(validate_orbital_backend(pa + pb, pa - pb, a, b, 4, 1, true).backend ==
              ElectronicBackend::AoP,
          "more occupied alpha MOs than available falls back to AO-P");
    check(validate_orbital_backend(pa + pb, pa - pb, a, b, 1, 1, false).backend ==
              ElectronicBackend::AoP,
          "unrestricted data treated as restricted falls back to AO-P");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)validate_orbital_backend(pa + pb, pa - pb, a, b, 1, 1, true, -1); },
        "a negative relative tolerance is rejected");
    return kumo_test::exit_code();
}
