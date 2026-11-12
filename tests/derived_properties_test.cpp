// Derived density properties of the Evaluator: ELF (total, alpha, beta), the
// reduced density gradient and sign(lambda2) rho.
//
//  - Each derived value equals its formula applied to the Evaluator's own
//    density, gradient, Hessian and tau. With AO_MO these are exactly the
//    inputs used; with AO_P the density values come from the grid kernel,
//    which agrees with the derivative kernel only to rounding, so the check
//    is restricted to rho >= 1e-4 and uses a looser bound.
//  - The Pauli kinetic energy density used by ELF agrees with the
//    cancellation-free orbital form (Lagrange identity)
//      tau_s - |grad rho_s|^2 / (8 rho_s)
//        = 1 / (2 rho_s) sum_{i<j} |phi_i grad phi_j - phi_j grad phi_i|^2.
//  - For a closed shell, the spin-resolved total ELF equals the
//    spin-unpolarized formula with the (3 pi^2)^(2/3) constant.
//  - Ranges, AO_P / AO_MO agreement, the values where there is no density,
//    invalid components, and an O-H bond critical point of water.
#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace chem::kumo;
using kumo_test::check;
using kumo_test::num;

constexpr double pi = 3.141592653589793238462643383279502884;

struct Point {
    double rho = 0.0;
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
    double tau = 0.0;
};

// Value, gradient, Hessian and tau of one component at point p.
Point at(const Results& r, DensityComponent c, std::size_t p) {
    const DensityDescriptor d{c};
    Point out;
    out.rho = r.get(d.value()).component()[p];
    for (std::size_t a = 0; a < 3; ++a)
        out.gradient[static_cast<Eigen::Index>(a)] = r.get(d.gradient()).component(a)[p];
    const auto& h = r.get(d.hessian());
    const std::size_t row[] = {0, 0, 0, 1, 1, 2}, col[] = {0, 1, 2, 1, 2, 2};
    for (std::size_t k = 0; k < 6; ++k)
        out.hessian(static_cast<Eigen::Index>(row[k]), static_cast<Eigen::Index>(col[k])) =
            out.hessian(static_cast<Eigen::Index>(col[k]), static_cast<Eigen::Index>(row[k])) =
                h.component(k)[p];
    out.tau = r.get(d.tau()).component()[p];
    return out;
}

double pauli(const Point& s) {
    return s.tau - s.gradient.squaredNorm() / (8.0 * s.rho);
}

double fermi(double rho) {
    return 0.3 * std::pow(6.0 * pi * pi, 2.0 / 3.0) * std::pow(rho, 5.0 / 3.0);
}

double elf_from(double d, double d0) {
    return 1.0 / (1.0 + (d / d0) * (d / d0));
}

// Pauli kinetic energy density of one spin channel from its occupied MOs.
double pauli_from_orbitals(const Molecule& m, const Eigen::Vector3d& point, Spin channel,
                           double& rho) {
    const auto n = static_cast<std::size_t>(m.occupied_count(channel));
    std::vector<OrbitalDerivatives> mo;
    for (std::size_t i = 1; i <= n; ++i)
        mo.push_back(m.orbital_derivatives(point, i, channel));
    rho = 0.0;
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        rho += mo[i].value * mo[i].value;
        for (std::size_t j = i + 1; j < n; ++j)
            sum += (mo[i].value * mo[j].gradient - mo[j].value * mo[i].gradient).squaredNorm();
    }
    return sum / (2.0 * rho);
}

Request full_request() {
    Request request;
    for (const DensityDescriptor& d : {density::total, density::alpha, density::beta})
        request.add(d.value()).add(d.gradient()).add(d.hessian()).add(d.tau()).add(d.elf());
    request.add(density::total.reduced_gradient()).add(density::total.sign_lambda2_density());
    return request;
}

void check_molecule(const std::string& name) {
    const auto molecule = Molecule::from_fchk(name);
    const auto grid = Grid::around(molecule.atoms(), 2.0, 0.45);
    const auto n = static_cast<std::size_t>(grid.counts.prod());
    const auto request = full_request();

    std::vector<Backend> backends{Backend::AO_P};
    if (molecule.electronic_structure_info().ao_mo_eligible()) backends.push_back(Backend::AO_MO);
    std::vector<Results> all;
    for (Backend backend : backends) {
        EvaluationOptions options;
        options.backend = backend;
        all.push_back(Evaluator(molecule, options).evaluate(grid, request));
    }

    for (std::size_t b = 0; b < all.size(); ++b) {
        const auto& r = all[b];
        const bool exact = backends[b] == Backend::AO_MO;
        const double floor = exact ? 0.0 : 1e-4;
        const std::string label = name + (exact ? " AO_MO" : " AO_P") + ": ";
        double formula_error = 0.0, closed_error = 0.0;
        bool ranges = true, sign_matches = true;
        for (std::size_t p = 0; p < n; ++p) {
            const auto t = at(r, DensityComponent::Total, p);
            const auto a = at(r, DensityComponent::Alpha, p);
            const auto be = at(r, DensityComponent::Beta, p);

            // The documented formulas applied to the Evaluator's own inputs.
            double d = 0.0, d0 = 0.0;
            for (const Point* s : {&a, &be})
                if (s->rho > 0.0) {
                    d += pauli(*s);
                    d0 += fermi(s->rho);
                }
            const double elf_total = r.get(density::total.elf()).component()[p];
            const double elf_alpha = r.get(density::alpha.elf()).component()[p];
            const double elf_beta = r.get(density::beta.elf()).component()[p];
            const double s = r.get(density::total.reduced_gradient()).component()[p];
            const double sl = r.get(density::total.sign_lambda2_density()).component()[p];
            const bool compare = t.rho > floor;
            if (compare && d0 > 0.0)
                formula_error = std::max(formula_error, std::abs(elf_total - elf_from(d, d0)));
            if (compare && a.rho > floor)
                formula_error =
                    std::max(formula_error, std::abs(elf_alpha - elf_from(pauli(a), fermi(a.rho))));
            if (compare && be.rho > floor)
                formula_error = std::max(formula_error,
                                         std::abs(elf_beta - elf_from(pauli(be), fermi(be.rho))));
            if (compare && t.rho > 0.0) {
                const double expected =
                    std::min(t.gradient.norm() /
                                 (2.0 * std::cbrt(3.0 * pi * pi) * std::pow(t.rho, 4.0 / 3.0)),
                             reduced_gradient_cap);
                formula_error = std::max(formula_error, std::abs(s - expected) / expected);
            }

            // Closed shell: spin-unpolarized ELF of the total density.
            if (!molecule.unrestricted() && t.rho > floor) {
                const double d_closed = t.tau - t.gradient.squaredNorm() / (8.0 * t.rho);
                const double d0_closed =
                    0.3 * std::pow(3.0 * pi * pi, 2.0 / 3.0) * std::pow(t.rho, 5.0 / 3.0);
                closed_error =
                    std::max(closed_error, std::abs(elf_total - elf_from(d_closed, d0_closed)));
            }

            ranges = ranges && elf_total >= 0.0 && elf_total <= 1.0 && elf_alpha >= 0.0 &&
                     elf_alpha <= 1.0 && elf_beta >= 0.0 && elf_beta <= 1.0 && s >= 0.0 &&
                     s <= reduced_gradient_cap &&
                     std::abs(std::abs(sl) - std::abs(t.rho)) <=
                         (exact ? 0.0 : 1e-12 * std::max(1.0, std::abs(t.rho)));

            // Sign of lambda2 from the closed-form (trigonometric) eigenvalues,
            // where it is not ambiguous.
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> direct;
            direct.computeDirect(t.hessian, Eigen::EigenvaluesOnly);
            const double lambda2 = direct.eigenvalues()[1];
            if (std::abs(lambda2) > 1e-8 * t.hessian.norm() && t.rho != 0.0)
                sign_matches = sign_matches && ((lambda2 < 0.0) == (sl * t.rho < 0.0));
        }
        check(formula_error < (exact ? 1e-12 : 1e-10),
              label + "derived values follow their formulas (max error " + num(formula_error) +
                  ")");
        if (!molecule.unrestricted())
            check(closed_error < 1e-10, label +
                                            "closed-shell ELF equals the spin-unpolarized "
                                            "formula (max error " +
                                            num(closed_error) + ")");
        check(ranges, label + "ELF in [0, 1], 0 <= s <= cap, |sign(lambda2) rho| == rho");
        check(sign_matches, label + "sign of sign(lambda2) rho matches lambda2");
    }

    // AO_P and AO_MO agree where the density is not tiny.
    if (all.size() == 2) {
        double difference = 0.0;
        for (std::size_t p = 0; p < n; ++p) {
            if (all[0].get(density::total.value()).component()[p] < 1e-4) continue;
            for (const Property& key :
                 {density::total.elf(), density::alpha.elf(), density::beta.elf()})
                difference = std::max(difference, std::abs(all[0].get(key).component()[p] -
                                                           all[1].get(key).component()[p]));
        }
        check(difference < 1e-8,
              name + ": AO_P and AO_MO ELF agree where rho >= 1e-4 (max " + num(difference) + ")");
    }

    // The total ELF alone (no Total derivatives evaluated, first AO
    // derivatives only) agrees with the one of the full request. It is not
    // bit-identical: AO screening also looks at the second-derivative panels
    // when the full request needs them.
    for (std::size_t b = 0; b < all.size(); ++b) {
        EvaluationOptions options;
        options.backend = backends[b];
        Request alone;
        alone.add(density::total.elf());
        const auto r = Evaluator(molecule, options).evaluate(grid, alone);
        const auto x = r.get(density::total.elf()).component();
        const auto y = all[b].get(density::total.elf()).component();
        double difference = 0.0;
        const auto rho = all[b].get(density::total.value()).component();
        for (std::size_t p = 0; p < n; ++p)
            if (rho[p] >= 1e-4) difference = std::max(difference, std::abs(x[p] - y[p]));
        check(difference < 1e-10,
              name +
                  ": a total-ELF-only request agrees with the full request where rho >= 1e-4 "
                  "(max " +
                  num(difference) + ")");
    }

    // Pauli kinetic energy density against its orbital form.
    double lagrange_error = 0.0;
    std::size_t lagrange_points = 0;
    const auto& r = all.front();
    for (std::size_t p = 0; p < n; p += 23) {
        const auto i =
            static_cast<int>(p / static_cast<std::size_t>(grid.counts[1] * grid.counts[2]));
        const auto j = static_cast<int>(p / static_cast<std::size_t>(grid.counts[2]) %
                                        static_cast<std::size_t>(grid.counts[1]));
        const auto k = static_cast<int>(p % static_cast<std::size_t>(grid.counts[2]));
        const Eigen::Vector3d point = grid.origin + grid.steps * Eigen::Vector3d(i, j, k);
        for (Spin channel : {Spin::Alpha, Spin::Beta}) {
            const auto c =
                channel == Spin::Alpha ? DensityComponent::Alpha : DensityComponent::Beta;
            const auto s = at(r, c, p);
            if (s.rho < 1e-3) continue;
            double rho = 0.0;
            const double orbital = pauli_from_orbitals(molecule, point, channel, rho);
            lagrange_error =
                std::max(lagrange_error, std::abs(pauli(s) - orbital) / std::max(s.tau, 1.0));
            ++lagrange_points;
        }
    }
    check(lagrange_points > 50, name + ": enough points for the orbital check");
    check(lagrange_error < 1e-9, name +
                                     ": Pauli kinetic energy density matches the orbital form "
                                     "(max relative error " +
                                     num(lagrange_error) + ")");

    // Far from the molecule every AO is cut off: no density, ELF 0, s capped.
    Grid far = grid;
    far.origin = grid.origin + Eigen::Vector3d(500.0, 0.0, 0.0);
    far.counts = Eigen::Vector3i(2, 2, 2);
    const auto empty = Evaluator(molecule).evaluate(far, request);
    bool defaults = true;
    for (std::size_t p = 0; p < 8; ++p)
        defaults =
            defaults && empty.get(density::total.value()).component()[p] == 0.0 &&
            empty.get(density::total.elf()).component()[p] == 0.0 &&
            empty.get(density::alpha.elf()).component()[p] == 0.0 &&
            empty.get(density::total.reduced_gradient()).component()[p] == reduced_gradient_cap &&
            empty.get(density::total.sign_lambda2_density()).component()[p] == 0.0;
    check(defaults, name + ": without density ELF is 0, s is the cap and sign(lambda2) rho is 0");

    // Components the derived properties are not defined for.
    const Evaluator evaluator(molecule);
    for (const Property& key :
         {density::spin.elf(), density::alpha.reduced_gradient(), density::spin.reduced_gradient(),
          density::beta.sign_lambda2_density(), density::spin.sign_lambda2_density()}) {
        Request bad;
        bad.add(key);
        kumo_test::check_throws<std::invalid_argument>(
            [&] { (void)evaluator.evaluate(grid, bad); },
            name + ": an undefined derived property is rejected");
    }
}

// Water: along an O-H bond, the reduced gradient vanishes at the bond
// critical point, where sign(lambda2) rho is negative (a bond) and rho is
// about 0.37 e/bohr^3 (B3LYP/def2-TZVP).
void check_water_bond(const std::string& name) {
    const auto molecule = Molecule::from_fchk(name);
    const auto& atoms = molecule.atoms();
    const auto oxygen = std::find_if(atoms.begin(), atoms.end(),
                                     [](const Atom& a) { return a.atomic_number == 8; });
    const auto hydrogen = std::find_if(atoms.begin(), atoms.end(),
                                       [](const Atom& a) { return a.atomic_number == 1; });
    if (!check(oxygen != atoms.end() && hydrogen != atoms.end(), name + ": water has O and H"))
        return;
    Grid line;
    const int count = 2001;
    line.origin = oxygen->position;
    line.steps = Eigen::Matrix3d::Zero();
    line.steps.col(0) = (hydrogen->position - oxygen->position) / (count - 1);
    line.steps(1, 1) = line.steps(2, 2) = 1.0;
    line.counts = Eigen::Vector3i(count, 1, 1);
    Request request;
    request.add(density::total.value())
        .add(density::total.reduced_gradient())
        .add(density::total.sign_lambda2_density())
        .add(density::total.elf());
    const auto r = Evaluator(molecule).evaluate(line, request);
    const auto s = r.get(density::total.reduced_gradient()).component();
    // Interior minimum (excluding the nuclei, where s is also small).
    std::size_t best = count / 10;
    for (std::size_t p = count / 10; p < static_cast<std::size_t>(count) * 9 / 10; ++p)
        if (s[p] < s[best]) best = p;
    const double rho = r.get(density::total.value()).component()[best];
    const double sl = r.get(density::total.sign_lambda2_density()).component()[best];
    const double elf = r.get(density::total.elf()).component()[best];
    std::cout << name
              << ": O-H bond critical point at t = " << static_cast<double>(best) / (count - 1)
              << ": s = " << s[best] << ", sign(lambda2) rho = " << sl << ", ELF = " << elf << '\n';
    check(s[best] < 0.02,
          name + ": s nearly vanishes at the O-H bond critical point (got " + num(s[best]) + ")");
    check(sl < 0.0 && rho > 0.33 && rho < 0.41,
          name + ": sign(lambda2) rho is about -0.37 at the O-H bond critical point (got " +
              num(sl) + ")");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 2, kumo_test::unlimited, "WATER.fchk FILE.fchk..."))
        return 2;
    return kumo_test::run([&] {
        check_water_bond(argv[1]);
        for (int i = 1; i < argc; ++i)
            check_molecule(argv[i]);
    });
}
