// Evaluator::evaluate at arbitrary points, and the orbital-energy fields
// (ALIE, LEA, LEAE, local electronegativity and hardness):
//  - every property at points equals its grid evaluation at the same points
//    (1e-12; 2e-8 for the electrostatic potential, evaluated in other
//    clusters), with both backends;
//  - the orbital-energy fields equal their formulas evaluated from
//    Molecule::orbital_derivatives() and the orbital energies, and the ALIE
//    equals average_local_ionization_energy();
//  - LEAE vanishes for a molecule without negative-energy virtual orbitals,
//    and the documented values where there is no density;
//  - identical results with 1, 2 and 4 threads, no points, invalid input.
//
// Usage: kumo-point-evaluation-test WATER.fchk O2.fchk ACROLEIN.fchk
#include "chem/kumo.hpp"
#include "omp_support.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace chem::kumo;
using kumo_test::check;
using kumo_test::num;
using Points = Eigen::Matrix<double, 3, Eigen::Dynamic>;

// The points of `grid` in Field order, computed as the library does.
Points points_of(const Grid& grid) {
    const int nx = grid.counts[0], ny = grid.counts[1], nz = grid.counts[2];
    Points points(3, static_cast<Eigen::Index>(nx) * ny * nz);
    Eigen::Index p = 0;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                points.col(p++) =
                    grid.origin + grid.steps * Eigen::Vector3d(static_cast<double>(i),
                                                               static_cast<double>(j),
                                                               static_cast<double>(k));
    return points;
}

// Every property defined for the molecule.
Request full_request(const Molecule& molecule) {
    Request request;
    const DensityDescriptor components[] = {density::total, density::alpha, density::beta,
                                            density::spin};
    for (const auto& d : components) {
        request.add(d.value()).add(d.gradient()).add(d.hessian()).add(d.laplacian()).add(d.tau());
        request.add(d.hamiltonian_kinetic_energy());
        if (d.component != DensityComponent::Spin) request.add(d.elf()).add(d.lol());
    }
    const auto& t = density::total;
    request.add(t.reduced_gradient())
        .add(t.sign_lambda2_density())
        .add(t.interaction_region_indicator())
        .add(t.density_overlap_regions_indicator())
        .add(t.potential_energy_density())
        .add(t.energy_density())
        .add(t.shannon_entropy_density())
        .add(t.fisher_information_density())
        .add(t.electrostatic_potential())
        .add(t.average_local_ionization_energy())
        .add(t.local_electron_affinity())
        .add(t.local_electron_attachment_energy())
        .add(t.local_electronegativity())
        .add(t.local_hardness())
        .add(density::spin.spin_polarization());
    const auto homo = static_cast<std::size_t>(molecule.occupied_count(Spin::Alpha));
    request.add(orbital::alpha(1).value()).add(orbital::alpha(homo).value());
    request.add(orbital::alpha(homo + 1).value());
    if (molecule.unrestricted())
        request.add(
            orbital::beta(static_cast<std::size_t>(molecule.occupied_count(Spin::Beta))).value());
    return request;
}

bool same_results(const Results& a, const Results& b, const Request& request) {
    for (const auto& key : request.properties())
        for (std::size_t c = 0; c < a.get(key).components.size(); ++c)
            if (!std::ranges::equal(a.get(key).component(c), b.get(key).component(c))) return false;
    return true;
}

void check_points_against_grid(const Molecule& molecule, const std::string& name) {
    const Grid grid = Grid::around(molecule.atoms(), 1.5, 0.7);
    const Points points = points_of(grid);
    const Request request = full_request(molecule);
    for (const Backend backend : {Backend::AO_P, Backend::AO_MO}) {
        const auto label = name + (backend == Backend::AO_P ? " AO-P" : " AO-MO");
        const Evaluator evaluator(molecule, {.backend = backend});
        const Results on_grid = evaluator.evaluate(grid, request);
        const Results at_points = evaluator.evaluate(points, request);
        check(at_points.backend_used() == on_grid.backend_used(),
              label + ": the same backend at points and on the grid");
        double worst = 0, worst_potential = 0;
        for (const auto& key : request.properties()) {
            const auto& expected = on_grid.get(key);
            const auto& actual = at_points.get(key);
            if (!check(actual.points == expected.points &&
                           actual.components.size() == expected.components.size(),
                       label + ": one value per point and component"))
                return;
            for (std::size_t c = 0; c < expected.components.size(); ++c)
                for (std::size_t p = 0; p < expected.points; ++p) {
                    const double e = expected.component(c)[p], a = actual.component(c)[p];
                    if (key.kind == PropertyKind::ElectrostaticPotential)
                        worst_potential = std::max(worst_potential, std::abs(a - e));
                    else
                        worst = std::max(worst, std::abs(a - e) / std::max(1.0, std::abs(e)));
                }
        }
        check(worst <= 1e-12,
              label + ": values at points equal grid values to 1e-12 (max " + num(worst) + ")");
        check(worst_potential <= 2e-8,
              label +
                  ": the electrostatic potential at points agrees with the grid within "
                  "2e-8 (max " +
                  num(worst_potential) + ")");
    }

    // Identical results with 1, 2 and 4 threads (the clusters of the
    // potential and the point blocks do not depend on them).
    kumo_test::set_threads(1);
    const Results serial = Evaluator(molecule).evaluate(points, request);
    for (int threads : {2, 4}) {
        kumo_test::set_threads(threads);
        check(same_results(serial, Evaluator(molecule).evaluate(points, request), request),
              name + ": identical results at points with " + std::to_string(threads) + " threads");
    }
    kumo_test::set_threads(1);
}

// The orbital-energy fields from orbital_derivatives() and the energies.
void check_orbital_energy_formulas(const Molecule& molecule, const std::string& name,
                                   std::size_t minimum_attachment_points) {
    const Grid grid = Grid::around(molecule.atoms(), 2.0, 1.1);
    const Points all = points_of(grid);
    Points points(3, 0);
    for (Eigen::Index p = 0; p < all.cols(); p += std::max<Eigen::Index>(1, all.cols() / 100)) {
        points.conservativeResize(Eigen::NoChange, points.cols() + 1);
        points.col(points.cols() - 1) = all.col(p);
    }
    const auto& t = density::total;
    Request request;
    request.add(t.average_local_ionization_energy())
        .add(t.local_electron_affinity())
        .add(t.local_electron_attachment_energy())
        .add(t.local_electronegativity())
        .add(t.local_hardness());
    const Results results = Evaluator(molecule).evaluate(points, request);
    const std::vector<double> alie_surface = average_local_ionization_energy(molecule, points);

    const bool unrestricted = molecule.unrestricted();
    const std::vector<Spin> channels =
        unrestricted ? std::vector<Spin>{Spin::Alpha, Spin::Beta} : std::vector<Spin>{Spin::Alpha};
    double worst = 0, worst_attachment = 0;
    std::size_t attachment_points = 0; // with |LEAE| >= 1e-3 hartree
    bool identical_alie = true;
    for (Eigen::Index p = 0; p < points.cols(); ++p) {
        double occupied = 0, occupied_energy = 0, virtual_density = 0, virtual_energy = 0,
               attachment = 0;
        for (const Spin channel : channels) {
            const auto& energies = molecule.orbital_energies(channel);
            const auto occupied_count = static_cast<std::size_t>(molecule.occupied_count(channel));
            for (std::size_t i = 0; i < molecule.orbital_count(channel); ++i) {
                const double phi =
                    molecule.orbital_derivatives(points.col(p), i + 1, channel).value;
                const double n = unrestricted ? 1.0 : 2.0;
                if (i < occupied_count) {
                    occupied += n * phi * phi;
                    occupied_energy += n * std::abs(energies[i]) * phi * phi;
                } else {
                    virtual_density += phi * phi;
                    virtual_energy += energies[i] * phi * phi;
                    if (energies[i] < 0) attachment += n * energies[i] * phi * phi;
                }
            }
        }
        const double alie = occupied_energy / occupied, lea = -virtual_energy / virtual_density;
        const double expected[] = {alie, lea, attachment / occupied, 0.5 * (alie + lea),
                                   0.5 * (alie - lea)};
        // LEAE (small, about 1e-2 hartree at most) also on its own scale.
        if (std::abs(expected[2]) >= 1e-3) {
            ++attachment_points;
            worst_attachment = std::max(worst_attachment,
                                        std::abs(results.get(t.local_electron_attachment_energy())
                                                     .component()[static_cast<std::size_t>(p)] -
                                                 expected[2]) /
                                            std::abs(expected[2]));
        }
        std::size_t f = 0;
        for (const auto& key : request.properties()) {
            const double actual = results.get(key).component()[static_cast<std::size_t>(p)];
            worst = std::max(worst, std::abs(actual - expected[f]) / std::max(1.0, std::abs(alie)));
            ++f;
        }
        identical_alie = identical_alie && results.get(t.average_local_ionization_energy())
                                                   .component()[static_cast<std::size_t>(p)] ==
                                               alie_surface[static_cast<std::size_t>(p)];
    }
    check(worst <= 1e-12, name +
                              ": ALIE, LEA, LEAE, chi and eta equal their orbital formulas to "
                              "1e-12 of the ALIE (max " +
                              num(worst) + ")");
    check(identical_alie, name + ": the ALIE equals average_local_ionization_energy()");
    // Molecules with negative-energy virtual orbitals (O2, acrolein): LEAE is
    // compared where it is substantial, relative to its own value.
    if (minimum_attachment_points > 0) {
        check(attachment_points >= minimum_attachment_points,
              name + ": at least " + std::to_string(minimum_attachment_points) +
                  " sampled points with |LEAE| >= 1e-3 hartree (got " +
                  std::to_string(attachment_points) + ")");
        check(worst_attachment <= 1e-12,
              name +
                  ": where |LEAE| >= 1e-3, it equals its orbital formula to 1e-12 relative "
                  "(max " +
                  num(worst_attachment) + ")");
    }

    // No density far from the molecule: every field is 0.
    Points far(3, 2);
    far.col(0) = Eigen::Vector3d(500, 0, 0);
    far.col(1) = Eigen::Vector3d(0, -800, 300);
    const Results empty = Evaluator(molecule).evaluate(far, request);
    bool zero = true;
    for (const auto& key : request.properties())
        for (const double v : empty.get(key).component())
            zero = zero && v == 0.0;
    check(zero, name + ": the orbital-energy fields are 0 where there is no density");
}

void check_input(const Molecule& water) {
    const Evaluator evaluator(water);
    Request request;
    request.add(density::total.value())
        .add(density::total.electrostatic_potential())
        .add(density::total.local_electron_attachment_energy())
        .add(orbital::alpha(1).value());
    const Results none = evaluator.evaluate(Points(3, 0), request);
    bool empty = true;
    for (const auto& key : request.properties())
        empty = empty && none.get(key).points == 0 && none.get(key).component().empty();
    check(empty, "no points give empty fields");

    Points bad(3, 2);
    bad.col(0) = Eigen::Vector3d(0, 0, 0);
    bad.col(1) = Eigen::Vector3d(std::numeric_limits<double>::quiet_NaN(), 0, 0);
    kumo_test::check_throws<std::invalid_argument>([&] { (void)evaluator.evaluate(bad, request); },
                                                   "a NaN point is rejected");
    Request orbital_request;
    orbital_request.add(orbital::alpha(10000).value());
    kumo_test::check_throws<std::out_of_range>(
        [&] { (void)evaluator.evaluate(Points(3, 1), orbital_request); },
        "a missing orbital is rejected at points");
    for (const Property key :
         {density::spin.average_local_ionization_energy(), density::alpha.local_electron_affinity(),
          density::beta.local_electron_attachment_energy(), density::spin.local_electronegativity(),
          density::alpha.local_hardness()}) {
        Request single;
        single.add(key);
        kumo_test::check_throws<std::invalid_argument>(
            [&] { (void)evaluator.evaluate(Points(3, 1), single); },
            "an orbital-energy field of a spin component is rejected");
    }

    // Water has no virtual orbital with negative energy: LEAE is 0.
    const Grid grid = Grid::around(water.atoms(), 1.0, 0.8);
    Request attachment;
    attachment.add(density::total.local_electron_attachment_energy())
        .add(density::total.local_electron_affinity());
    const Results r = evaluator.evaluate(grid, attachment);
    const auto leae = r.get(density::total.local_electron_attachment_energy()).component();
    const auto lea = r.get(density::total.local_electron_affinity()).component();
    check(std::all_of(leae.begin(), leae.end(), [](double v) { return v == 0.0; }),
          "water: LEAE is 0 everywhere (no virtual orbital with negative energy)");
    check(std::all_of(lea.begin(), lea.end(), [](double v) { return v < 0.0; }),
          "water: LEA is negative everywhere on the grid (virtual energies are positive)");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 3, 3, "WATER.fchk O2.fchk ACROLEIN.fchk")) return 2;
    return kumo_test::run([&] {
        const char* names[] = {"water", "O2", "acrolein"};
        for (int m = 0; m < 3; ++m) {
            const auto molecule = Molecule::from_fchk(argv[m + 1]);
            check_points_against_grid(molecule, names[m]);
            // Water has no virtual orbital with negative energy.
            check_orbital_energy_formulas(molecule, names[m], m == 0 ? 0 : 20);
            if (m == 0) check_input(molecule);
        }
    });
}
