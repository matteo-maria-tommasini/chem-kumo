// Electrostatic potential (detail::ElectrostaticPotential) against
// independent references:
//  - closed forms for s-function densities on one and on two centres;
//  - the multipole expansion far from the molecule, with the charge, dipole
//    and quadrupole of the electron density from the analytic AO moments
//    (BasisEvaluator::cartesian_moments, a separate code path);
//  - the Poisson equation, laplacian(V_el) = 4 pi rho, by finite differences,
//    with rho from the density evaluator; this covers every shell type;
//  - the nuclear term, and the potential at a nucleus (its own term left out).
#include "chem/kumo.hpp"
#include "chem/kumo/basis_evaluator.hpp"
#include "electrostatic_potential.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace chem::kumo;
using detail::ElectrostaticPotential;

constexpr double pi = std::numbers::pi;

// Potential of the product of two normalized contracted s functions with the
// given exponents and coefficients, centred at A and B, at `point`:
//   sum_pq w_p w_q K_pq (pi / p)^(3/2) erf(sqrt(p) |r - P|) / |r - P|.
double s_pair_potential(const Shell& a, const Shell& b, const Eigen::Vector3d& point) {
    double sum = 0.0;
    for (std::size_t p = 0; p < a.exponents.size(); ++p)
        for (std::size_t q = 0; q < b.exponents.size(); ++q) {
            const double alpha = a.exponents[p], beta = b.exponents[q], s = alpha + beta;
            const double wa = a.coefficients[p] * std::pow(2.0 * alpha / pi, 0.75);
            const double wb = b.coefficients[q] * std::pow(2.0 * beta / pi, 0.75);
            const Eigen::Vector3d center = (alpha * a.origin + beta * b.origin) / s;
            const double K = std::exp(-alpha * beta / s * (a.origin - b.origin).squaredNorm());
            const double r = (point - center).norm();
            sum += wa * wb * K * std::pow(pi / s, 1.5) * std::erf(std::sqrt(s) * r) / r;
        }
    return sum;
}

void closed_forms(const std::string& label, const BasisSet& basis) {
    // First AO of each shell, in FCHK order.
    std::vector<Eigen::Index> first;
    Eigen::Index offset = 0;
    for (const auto& shell : basis.shells()) {
        first.push_back(offset);
        offset += static_cast<Eigen::Index>(shell.function_count());
    }
    const auto& shells = basis.shells();
    const auto is_s = [](const Shell& s) { return s.angular_momentum == 0 || s.is_sp; };
    std::size_t other = shells.size();
    for (std::size_t k = 1; k < shells.size(); ++k)
        if (is_s(shells[k]) && shells[k].center != shells[0].center) {
            other = k;
            break;
        }
    if (!kumo_test::check(is_s(shells[0]) && other < shells.size(),
                          label + ": the basis has s shells on two atoms"))
        return;

    const auto n = static_cast<Eigen::Index>(basis.function_count());
    const Eigen::Index i = first[0], j = first[other];
    Eigen::MatrixXd diagonal = Eigen::MatrixXd::Zero(n, n);
    diagonal(i, i) = 1.0;
    Eigen::MatrixXd cross = Eigen::MatrixXd::Zero(n, n);
    cross(i, j) = cross(j, i) = 1.0;
    const ElectrostaticPotential one(basis, diagonal, {});
    const ElectrostaticPotential two(basis, cross, {});

    double worst = 0.0;
    const Eigen::Vector3d points[] = {
        shells[0].origin + Eigen::Vector3d(0.3, -0.2, 0.1),
        shells[0].origin + Eigen::Vector3d(1e-9, 0.0, 0.0),
        0.5 * (shells[0].origin + shells[other].origin) + Eigen::Vector3d(0.0, 0.7, 0.0),
        shells[other].origin + Eigen::Vector3d(-1.5, 2.0, 0.4), Eigen::Vector3d(8.0, -6.0, 3.0)};
    for (const auto& point : points) {
        const double expected_one = -s_pair_potential(shells[0], shells[0], point);
        const double expected_two = -2.0 * s_pair_potential(shells[0], shells[other], point);
        worst = std::max({worst, std::abs(one.electronic(point) - expected_one),
                          std::abs(two.electronic(point) - expected_two)});
    }
    std::cout << label << ": s closed forms, max difference " << worst << '\n';
    kumo_test::check(worst < 1e-12, label + ": s-function densities match the erf closed form (" +
                                        kumo_test::num(worst) + ")");
}

void far_field(const std::string& label, const BasisSet& basis, const Molecule& molecule) {
    const Eigen::MatrixXd P = molecule.density_matrix(DensityComponent::Total);
    const auto moments = BasisEvaluator(basis).cartesian_moments();
    const double electrons = (P.array() * moments.overlap.array()).sum();
    Eigen::Vector3d dipole;
    for (int a = 0; a < 3; ++a)
        dipole[a] = (P.array() * moments.first[static_cast<std::size_t>(a)].array()).sum();
    // second: xx, yy, zz, xy, xz, yz
    Eigen::Matrix3d second;
    const int map[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            second(a, b) =
                (P.array() * moments.second[static_cast<std::size_t>(map[a][b])].array()).sum();

    const ElectrostaticPotential potential(basis, P, {});
    double worst = 0.0, quadrupole = 0.0;
    const Eigen::Vector3d directions[] = {{1, 0, 0}, {0, -1, 0}, {0, 0, 1}, {0.6, 0.48, -0.64}};
    for (const auto& direction : directions) {
        const Eigen::Vector3d r = 1000.0 * direction.normalized();
        const double d = r.norm();
        // 1/|r - x| = 1/d + r.x / d^3 + (3 (r.x)^2 - d^2 |x|^2) / (2 d^5) + O(d^-4)
        const double quadrupole_term =
            (3.0 * r.dot(second * r) - d * d * second.trace()) / (2.0 * std::pow(d, 5));
        quadrupole = std::max(quadrupole, std::abs(quadrupole_term));
        const double expansion = electrons / d + r.dot(dipole) / (d * d * d) + quadrupole_term;
        worst = std::max(worst, std::abs(potential.electronic(r) + expansion));
    }
    // The expansion stops at the quadrupole; the residual is the octupole
    // and higher terms, (extent / d) times smaller or less. Requiring it to
    // be below 1% of the quadrupole term checks the charge and dipole far
    // more tightly and the quadrupole to 1%. Observed ratios: 2e-3 (water),
    // 1e-6 (O2), 2e-4 (acrolein), 7e-5 (HBC).
    const double tolerance = 1e-2 * quadrupole;
    std::cout << label << ": far field (1000 bohr), max difference from the multipole expansion "
              << worst << " (quadrupole term " << quadrupole << ")\n";
    kumo_test::check(worst < tolerance, label +
                                            ": far field matches the charge, dipole and quadrupole "
                                            "of the AO moments (" +
                                            kumo_test::num(worst) + ")");
}

void poisson(const std::string& label, const BasisSet& basis, const Molecule& molecule) {
    const ElectrostaticPotential potential(basis, molecule.density_matrix(DensityComponent::Total),
                                           {});
    const auto& atoms = molecule.atoms();
    // Points between and around the first atoms (at most 8, which bounds the
    // cost for large molecules), away from the nuclear cusps.
    std::vector<Eigen::Vector3d> points;
    for (std::size_t a = 0; a < std::min<std::size_t>(atoms.size(), 8); ++a) {
        points.push_back(atoms[a].position + Eigen::Vector3d(0.35, -0.25, 0.3));
        points.push_back(atoms[a].position + Eigen::Vector3d(-0.9, 0.6, 1.1));
        if (a + 1 < atoms.size())
            points.push_back(0.5 * (atoms[a].position + atoms[a + 1].position) +
                             Eigen::Vector3d(0.05, 0.1, -0.07));
    }
    // The difference is the O(h^4) truncation error of the stencil: it was
    // 1.2e-6 with h = 0.01 and 7.6e-8 with h = 0.005 (16 times smaller).
    const double h = 0.005;
    double worst = 0.0;
    for (const auto& r : points) {
        // Fourth-order central differences of the second derivatives.
        double laplacian = 0.0;
        const double centre = potential.electronic(r);
        for (int a = 0; a < 3; ++a) {
            Eigen::Vector3d e = Eigen::Vector3d::Zero();
            e[a] = h;
            laplacian += (-potential.electronic(r + 2 * e) + 16.0 * potential.electronic(r + e) -
                          30.0 * centre + 16.0 * potential.electronic(r - e) -
                          potential.electronic(r - 2 * e)) /
                         (12.0 * h * h);
        }
        const double expected = 4.0 * pi * molecule.density_derivatives(r).value;
        worst = std::max(worst, std::abs(laplacian - expected) / std::max(1.0, expected));
    }
    std::cout << label << ": Poisson equation at " << points.size()
              << " points, max relative difference " << worst << '\n';
    kumo_test::check(worst < 5e-7, label + ": laplacian(V_el) = 4 pi rho by finite differences (" +
                                       kumo_test::num(worst) + ")");
}

void nuclear_term(const std::string& label, const BasisSet& basis, const Molecule& molecule) {
    const ElectrostaticPotential potential(basis, molecule.density_matrix(DensityComponent::Total),
                                           molecule.atoms());
    const Eigen::Vector3d point(0.4, -1.3, 2.2);
    double expected = 0.0;
    for (const auto& atom : molecule.atoms())
        expected += atom.nuclear_charge / (atom.position - point).norm();
    kumo_test::check(std::abs(potential.nuclear(point) - expected) < 1e-14 * expected,
                     label + ": the nuclear term is sum Z / |r - R|");
    kumo_test::check(
        std::abs(potential(point) - potential.nuclear(point) - potential.electronic(point)) == 0.0,
        label + ": the total is nuclear + electronic");

    const Eigen::Vector3d nucleus = molecule.atoms().front().position;
    double others = 0.0;
    for (std::size_t a = 1; a < molecule.atoms().size(); ++a)
        others +=
            molecule.atoms()[a].nuclear_charge / (molecule.atoms()[a].position - nucleus).norm();
    kumo_test::check(potential.nuclear(nucleus) == others &&
                         std::abs(potential.nuclear(nucleus + Eigen::Vector3d(0.0, 5e-9, 0.0)) -
                                  others) < 1e-7 * others,
                     label + ": on a nucleus (within 1e-8 bohr) its own term is left out");
    kumo_test::check(std::isfinite(potential(nucleus)) && potential(nucleus) < 0.0,
                     label + ": the potential at a nucleus (others plus electrons) is finite and "
                             "negative");
    const double near = potential(nucleus + Eigen::Vector3d(0.0, 0.0, 1e-6));
    kumo_test::check(std::isfinite(near) && near > 1e5,
                     label + ": 1e-6 bohr from a nucleus the potential is large and finite");
    kumo_test::check(std::isfinite(potential.electronic(nucleus)),
                     label + ": the electronic term is finite at a nucleus");

    Eigen::Matrix<double, 3, Eigen::Dynamic> block(3, 2);
    block.col(0) = point;
    block.col(1) = point + Eigen::Vector3d(0.1, 0.2, 0.3);
    std::vector<double> values(2);
    potential.evaluate(block, values);
    kumo_test::check(values[0] == potential(block.col(0)) && values[1] == potential(block.col(1)),
                     label + ": block evaluation equals pointwise evaluation");
}

// The Evaluator returns the engine's values on every grid point, in Field
// order, for any block size (and, in OpenMP builds, any thread count), with
// or without other properties in the same request.
void evaluator(const std::string& label, const BasisSet& basis, const Molecule& molecule) {
    const ElectrostaticPotential potential(basis, molecule.density_matrix(DensityComponent::Total),
                                           molecule.atoms());
    Grid grid;
    grid.origin = Eigen::Vector3d(-3.1, -2.9, -3.3);
    grid.steps = Eigen::Matrix3d::Identity() * 0.61;
    grid.steps(0, 1) = 0.05;
    grid.counts = Eigen::Vector3i(9, 8, 11);
    const auto key = density::total.electrostatic_potential();

    Request request;
    request.add(key).add(density::total.value()).add(density::total.elf());
    std::size_t last_done = 0, last_total = 0;
    const Results results = Evaluator(molecule, {.block_size = 17,
                                                 .potential_tolerance = 0.0,
                                                 .progress =
                                                     [&](std::size_t done, std::size_t total) {
                                                         last_done = done;
                                                         last_total = total;
                                                     }})
                                .evaluate(grid, request);
    const auto values = results.get(key).component();
    bool same = values.size() == 9u * 8u * 11u;
    for (std::size_t p = 0; same && p < values.size(); ++p) {
        const std::size_t k = p % 11, j = (p / 11) % 8, i = p / 88;
        const Eigen::Vector3d point =
            grid.origin + grid.steps * Eigen::Vector3d(double(i), double(j), double(k));
        same = values[p] == potential(point);
    }
    kumo_test::check(same, label + ": with tolerance 0, Evaluator values equal the exact engine "
                                   "at every grid point");
    kumo_test::check(last_total > 0 && last_done == last_total,
                     label + ": progress reaches the total");

    Request alone;
    alone.add(key);
    const Results single =
        Evaluator(molecule, {.block_size = 1000, .potential_tolerance = 0.0}).evaluate(grid, alone);
    const auto again = single.get(key).component();
    kumo_test::check(std::equal(values.begin(), values.end(), again.begin(), again.end()),
                     label + ": the result does not depend on the block size or the request");

    // Default tolerance: within 1e-8 of the exact values, independently of the
    // block size and of the other properties requested.
    const Results fast = Evaluator(molecule).evaluate(grid, request);
    const Results fast_alone = Evaluator(molecule, {.block_size = 5}).evaluate(grid, alone);
    const auto approximate = fast.get(key).component();
    const auto approximate_alone = fast_alone.get(key).component();
    double worst = 0.0;
    for (std::size_t p = 0; p < values.size(); ++p)
        worst = std::max(worst, std::abs(approximate[p] - values[p]));
    std::cout << label << ": Evaluator with the default tolerance, max difference " << worst
              << '\n';
    kumo_test::check(worst <= 1e-8,
                     label + ": with the default tolerance, within 1e-8 of the exact values (" +
                         kumo_test::num(worst) + ")");
    kumo_test::check(std::equal(approximate.begin(), approximate.end(), approximate_alone.begin(),
                                approximate_alone.end()),
                     label + ": the approximate result does not depend on the block size or the "
                             "request");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { Evaluator(molecule, {.potential_tolerance = -1e-9}); },
        label + ": a negative potential tolerance is rejected");

    for (const auto component : {density::alpha, density::beta, density::difference}) {
        Request wrong;
        wrong.add(component.electrostatic_potential());
        kumo_test::check_throws<std::invalid_argument>(
            [&] { (void)Evaluator(molecule).evaluate(grid, wrong); },
            label + ": the electrostatic potential of a spin component is rejected");
    }
    Grid on_nucleus;
    on_nucleus.origin = molecule.atoms().back().position - Eigen::Vector3d(0.5, 0.5, 0.5);
    on_nucleus.steps = Eigen::Matrix3d::Identity() * 0.5;
    on_nucleus.counts = Eigen::Vector3i(3, 3, 3);
    const Results through = Evaluator(molecule).evaluate(on_nucleus, alone);
    kumo_test::check(std::abs(through.get(key).component()[13] -
                              potential(molecule.atoms().back().position)) <= 1e-8,
                     label + ": on a grid through a nucleus, the point on it gets the potential "
                             "at the nucleus");
}

// Far-field plans and cluster evaluation: the error is within the tolerance,
// some pairs are expanded, tolerance 0 is exact, and invalid use is rejected.
void clusters(const std::string& label, const BasisSet& basis, const Molecule& molecule) {
    const ElectrostaticPotential potential(basis, molecule.density_matrix(DensityComponent::Total),
                                           molecule.atoms());
    const double radius = 1.2;
    // Points on a sphere and inside it, around centres in and far from the molecule.
    const auto& atoms = molecule.atoms();
    const Eigen::Vector3d centres[] = {atoms.front().position + Eigen::Vector3d(0.3, 0.2, -0.4),
                                       atoms.back().position + Eigen::Vector3d(-2.0, 1.0, 1.5),
                                       Eigen::Vector3d(9.0, -7.0, 4.0)};
    Eigen::Matrix<double, 3, Eigen::Dynamic> offsets(3, 27);
    int q = 0;
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = -1; k <= 1; ++k, ++q) {
                // The centre, then directions alternately on the sphere and at 0.6 of it.
                const double scale = (i == 0 && j == 0 && k == 0) ? 0.0 : (q % 2 ? 1.0 : 0.6);
                offsets.col(q) = Eigen::Vector3d(i, j, k).normalized() * radius * scale;
            }
    std::vector<double> exact(27), approximate(27);
    for (const double tolerance : {0.0, 1e-4, 1e-8}) {
        const auto plan = potential.far_field_plan(radius, tolerance);
        const auto expanded =
            std::count_if(plan.level_distance2.begin(), plan.level_distance2.end(),
                          [](float d2) { return std::isfinite(d2); }) +
            static_cast<long>(plan.clusters.size());
        kumo_test::check(tolerance == 0.0 ? expanded == 0 : expanded > 0,
                         label + ": tolerance " + kumo_test::num(tolerance) +
                             (tolerance == 0.0 ? " expands no pair" : " expands some pairs"));
        double worst = 0.0;
        for (const auto& centre : centres) {
            Eigen::Matrix<double, 3, Eigen::Dynamic> points = offsets.colwise() + centre;
            potential.evaluate_cluster(points, centre, plan, approximate);
            for (Eigen::Index p = 0; p < points.cols(); ++p) {
                exact[static_cast<std::size_t>(p)] = potential(points.col(p));
                worst = std::max(worst, std::abs(approximate[static_cast<std::size_t>(p)] -
                                                 exact[static_cast<std::size_t>(p)]));
            }
        }
        std::cout << label << ": clusters with tolerance " << tolerance << ", max difference "
                  << worst << '\n';
        kumo_test::check(worst <= tolerance, label + ": cluster values within the tolerance " +
                                                 kumo_test::num(tolerance) + " (" +
                                                 kumo_test::num(worst) + ")");
    }
    const auto plan = potential.far_field_plan(radius, 1e-8);
    Eigen::Matrix<double, 3, Eigen::Dynamic> outside(3, 1);
    outside.col(0) = centres[0] + Eigen::Vector3d(0.0, 0.0, 1.3);
    std::vector<double> one(1);
    kumo_test::check_throws<std::invalid_argument>(
        [&] { potential.evaluate_cluster(outside, centres[0], plan, one); },
        label + ": a point outside the planned radius is rejected");
    ElectrostaticPotential::FarFieldPlan wrong = plan;
    wrong.level_distance2.pop_back();
    kumo_test::check_throws<std::invalid_argument>(
        [&] {
            potential.evaluate_cluster(offsets.colwise() + centres[0], centres[0], wrong, exact);
        },
        label + ": a plan for another potential is rejected");
    kumo_test::check_throws<std::invalid_argument>(
        [&] { (void)potential.far_field_plan(radius, -1.0); },
        label + ": a negative tolerance is rejected");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        for (int k = 1; k < argc; ++k) {
            const std::string label = std::filesystem::path(argv[k]).filename().string();
            const FormattedCheckpoint fchk(argv[k]);
            const BasisSet basis(fchk);
            const auto molecule = Molecule::from_fchk(argv[k]);
            closed_forms(label, basis);
            far_field(label, basis, molecule);
            poisson(label, basis, molecule);
            nuclear_term(label, basis, molecule);
            evaluator(label, basis, molecule);
            clusters(label, basis, molecule);
            const auto n = static_cast<Eigen::Index>(basis.function_count());
            kumo_test::check_throws<std::invalid_argument>(
                [&] { ElectrostaticPotential(basis, Eigen::MatrixXd::Zero(n + 1, n + 1), {}); },
                label + ": a density matrix of the wrong size is rejected");
            Eigen::MatrixXd bad = Eigen::MatrixXd::Zero(n, n);
            bad(0, 0) = std::nan("");
            kumo_test::check_throws<std::invalid_argument>(
                [&] { ElectrostaticPotential(basis, bad, {}); },
                label + ": a nonfinite density matrix is rejected");
        }
    });
}
