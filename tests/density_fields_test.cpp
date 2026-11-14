// Density fields of the Evaluator added in 1.12.0: LOL, IRI, DORI, the
// kinetic and energy densities K, V and H, the spin polarization, and the
// Shannon and Fisher densities.
//
//  - Each value equals its formula applied to the Evaluator's own density,
//    gradient, Hessian, Laplacian and tau of the same request.
//  - K agrees with the orbital form -(1/2) sum_i phi_i lap phi_i, evaluated
//    from Molecule::orbital_derivatives(), and DORI with central finite
//    differences of |grad rho|^2 / rho^2 from Molecule::density_derivatives().
//  - For a closed shell the spin-resolved total LOL equals the
//    spin-unpolarized formula with the (3 pi^2)^(2/3) constant, and the spin
//    polarization vanishes.
//  - Ranges, the values where there is no density, invalid components, and
//    identical results with 1, 2 and 4 threads.
#include "chem/kumo.hpp"
#include "omp_support.hpp"
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
using kumo_test::check;
using kumo_test::num;

constexpr double pi = std::numbers::pi;

const DensityDescriptor all_components[] = {density::total, density::alpha, density::beta,
                                            density::difference};

Request full_request() {
    Request request;
    for (const auto& d : all_components)
        request.add(d.value())
            .add(d.gradient())
            .add(d.hessian())
            .add(d.laplacian())
            .add(d.tau())
            .add(d.hamiltonian_kinetic_energy());
    for (const auto& d : {density::total, density::alpha, density::beta})
        request.add(d.lol());
    request.add(density::total.interaction_region_indicator())
        .add(density::total.density_overlap_regions_indicator())
        .add(density::total.potential_energy_density())
        .add(density::total.energy_density())
        .add(density::total.shannon_entropy_density())
        .add(density::total.fisher_information_density())
        .add(density::difference.spin_polarization());
    return request;
}

double at(const Results& r, Property key, std::size_t p, std::size_t component = 0) {
    return r.get(key).component(component)[p];
}

Eigen::Vector3d gradient_at(const Results& r, const DensityDescriptor& d, std::size_t p) {
    return {at(r, d.gradient(), p, 0), at(r, d.gradient(), p, 1), at(r, d.gradient(), p, 2)};
}

// Evaluator Hessian order: xx, xy, xz, yy, yz, zz.
Eigen::Matrix3d hessian_at(const Results& r, const DensityDescriptor& d, std::size_t p) {
    const auto h = [&](std::size_t c) { return at(r, d.hessian(), p, c); };
    Eigen::Matrix3d m;
    m << h(0), h(1), h(2), h(1), h(3), h(4), h(2), h(4), h(5);
    return m;
}

// LOL from (rho_s, tau_s) of the channels given, with the spin-resolved
// Thomas-Fermi constant.
double lol_formula(std::initializer_list<std::pair<double, double>> channels) {
    const double fermi = 0.3 * std::pow(6.0 * pi * pi, 2.0 / 3.0);
    double tau = 0.0, uniform = 0.0;
    for (const auto& [rho, t] : channels) {
        tau += t;
        if (rho > 0.0) uniform += fermi * std::pow(rho, 5.0 / 3.0);
    }
    if (!(uniform > 0.0)) return 0.0;
    if (!(tau > 0.0)) return 1.0;
    const double ratio = uniform / tau;
    return ratio / (1.0 + ratio);
}

// Largest |a - b| / (scale + floor) over the points, for the report.
struct Worst {
    double ratio = 0.0;

    void add(double a, double b, double scale) { ratio = std::max(ratio, std::abs(a - b) / scale); }
};

void check_molecule(const std::string& name) {
    const auto molecule = Molecule::from_fchk(name);
    const Grid grid = Grid::around(molecule.atoms(), 1.5, 0.45);
    const std::size_t n = static_cast<std::size_t>(grid.counts.prod());
    const Request request = full_request();
    kumo_test::set_threads(1);
    const Results r = Evaluator(molecule).evaluate(grid, request);
    const bool closed_shell =
        !molecule.unrestricted() &&
        molecule.occupied_count(SpinChannel::Alpha) == molecule.occupied_count(SpinChannel::Beta);

    // Formulas applied to the Evaluator's own inputs, to rounding.
    Worst energy, lol, iri, dori, polarization, shannon, fisher;
    bool ranges = true;
    for (std::size_t p = 0; p < n; ++p) {
        for (const auto& d : all_components) {
            const double tau = at(r, d.tau(), p), lap = at(r, d.laplacian(), p);
            energy.add(at(r, d.hamiltonian_kinetic_energy(), p), tau - 0.25 * lap,
                       std::abs(tau) + 0.25 * std::abs(lap) + 1e-300);
        }
        const double tau = at(r, density::total.tau(), p);
        const double lap = at(r, density::total.laplacian(), p);
        const double scale = std::abs(tau) + 0.25 * std::abs(lap) + 1e-300;
        energy.add(at(r, density::total.potential_energy_density(), p), 0.25 * lap - 2.0 * tau,
                   scale);
        energy.add(at(r, density::total.energy_density(), p), 0.25 * lap - tau, scale);

        const double rho = at(r, density::total.value(), p);
        const double ra = at(r, density::alpha.value(), p), rb = at(r, density::beta.value(), p);
        const double ta = at(r, density::alpha.tau(), p), tb = at(r, density::beta.tau(), p);
        lol.add(at(r, density::total.lol(), p), lol_formula({{ra, ta}, {rb, tb}}), 1.0);
        lol.add(at(r, density::alpha.lol(), p), lol_formula({{ra, ta}}), 1.0);
        lol.add(at(r, density::beta.lol(), p), lol_formula({{rb, tb}}), 1.0);

        const Eigen::Vector3d g = gradient_at(r, density::total, p);
        const double expected_iri =
            rho > 0.0 ? std::min(g.norm() / std::pow(rho, 1.1), interaction_region_indicator_cap)
                      : interaction_region_indicator_cap;
        iri.add(at(r, density::total.interaction_region_indicator(), p), expected_iri,
                std::abs(expected_iri) + 1e-300);

        if (rho > 0.0 && g.squaredNorm() > 0.0) {
            // theta = |grad f|^2 / f^3, f = |g|^2 / rho^2, grad f = 2 (H g rho - |g|^2 g) / rho^3.
            const Eigen::Matrix3d h = hessian_at(r, density::total, p);
            const double f = g.squaredNorm() / (rho * rho);
            const Eigen::Vector3d grad_f =
                2.0 * (h * g * rho - g.squaredNorm() * g) / (rho * rho * rho);
            const double theta = grad_f.squaredNorm() / (f * f * f);
            if (std::isfinite(theta))
                dori.add(at(r, density::total.density_overlap_regions_indicator(), p),
                         theta / (1.0 + theta), 1.0);
        }

        const double zeta = at(r, density::difference.spin_polarization(), p);
        polarization.add(zeta, rho > 0.0 ? at(r, density::difference.value(), p) / rho : 0.0, 1.0);
        ranges = ranges && std::abs(zeta) <= 1.0 + 1e-12 && (!closed_shell || zeta == 0.0);
        for (const auto& key : {density::total.lol(), density::alpha.lol(), density::beta.lol(),
                                density::total.density_overlap_regions_indicator()}) {
            const double v = at(r, key, p);
            ranges = ranges && v >= 0.0 && v <= 1.0;
        }
        shannon.add(at(r, density::total.shannon_entropy_density(), p),
                    rho > 0.0 ? -rho * std::log(rho) : 0.0,
                    rho > 0.0 ? rho * (1.0 + std::abs(std::log(rho))) : 1e-300);
        fisher.add(at(r, density::total.fisher_information_density(), p),
                   rho > 0.0 ? g.squaredNorm() / rho : 0.0,
                   rho > 0.0 ? g.squaredNorm() / rho + 1e-300 : 1e-300);
    }
    const auto formula = [&](const char* what, const Worst& w) {
        check(w.ratio <= 1e-13,
              name + ": " + what + " matches its formula to 1e-13 (worst " + num(w.ratio) + ")");
    };
    formula("K, V and H", energy);
    formula("LOL", lol);
    formula("IRI", iri);
    formula("the spin polarization", polarization);
    formula("the Shannon density", shannon);
    formula("the Fisher density", fisher);
    check(dori.ratio <= 1e-12,
          name + ": DORI matches its formula to 1e-12 (worst " + num(dori.ratio) + ")");
    check(ranges, name + ": LOL and DORI lie in [0, 1], |zeta| <= 1" +
                      (closed_shell ? ", and zeta = 0 for a closed shell" : std::string()));

    // Closed shell: the spin-resolved total LOL is the unpolarized formula.
    if (closed_shell) {
        const double fermi = 0.3 * std::pow(3.0 * pi * pi, 2.0 / 3.0);
        Worst unpolarized;
        for (std::size_t p = 0; p < n; ++p) {
            const double rho = at(r, density::total.value(), p);
            const double tau = at(r, density::total.tau(), p);
            if (!(rho > 0.0 && tau > 0.0)) continue;
            const double ratio = fermi * std::pow(rho, 5.0 / 3.0) / tau;
            unpolarized.add(at(r, density::total.lol(), p), ratio / (1.0 + ratio), 1.0);
        }
        check(unpolarized.ratio <= 1e-12, name +
                                              ": closed-shell LOL matches the unpolarized form "
                                              "(worst " +
                                              num(unpolarized.ratio) + ")");
    }

    // Independent checks at points spread over the grid with enough density.
    Worst orbital_k, finite_difference;
    std::size_t sampled = 0;
    const std::size_t stride = std::max<std::size_t>(1, n / 200);
    for (std::size_t p = 0; p < n; p += stride) {
        const auto ijk = Eigen::Vector3i(static_cast<int>(p / (grid.counts[1] * grid.counts[2])),
                                         static_cast<int>(p / grid.counts[2] % grid.counts[1]),
                                         static_cast<int>(p % grid.counts[2]));
        const Eigen::Vector3d x = grid.origin + grid.steps * ijk.cast<double>();
        const double rho = at(r, density::total.value(), p);
        if (rho < 1e-3) continue;
        ++sampled;
        // K_s = -(1/2) sum_i phi_i lap phi_i over the occupied orbitals of each channel.
        for (const auto& [channel, descriptor] : {std::pair{SpinChannel::Alpha, density::alpha},
                                                  std::pair{SpinChannel::Beta, density::beta}}) {
            double k = 0.0;
            for (std::size_t i = 1; i <= molecule.occupied_count(channel); ++i) {
                const auto phi = molecule.orbital_derivatives(x, i, channel);
                k -= 0.5 * phi.value * phi.hessian.trace();
            }
            const double tau = at(r, descriptor.tau(), p), lap = at(r, descriptor.laplacian(), p);
            orbital_k.add(at(r, descriptor.hamiltonian_kinetic_energy(), p), k,
                          std::abs(tau) + 0.25 * std::abs(lap));
        }
        // theta from central differences (step 1e-4 bohr) of f = |grad rho|^2 / rho^2.
        const auto f = [&](const Eigen::Vector3d& y) {
            const auto d = molecule.density_derivatives(y);
            return d.gradient.squaredNorm() / (d.value * d.value);
        };
        const double h = 1e-4, f0 = f(x);
        if (!(f0 > 1e-6)) continue;
        Eigen::Vector3d grad_f;
        for (int a = 0; a < 3; ++a) {
            Eigen::Vector3d e = Eigen::Vector3d::Zero();
            e[a] = h;
            grad_f[a] = (f(x + e) - f(x - e)) / (2.0 * h);
        }
        const double theta = grad_f.squaredNorm() / (f0 * f0 * f0);
        finite_difference.add(at(r, density::total.density_overlap_regions_indicator(), p),
                              theta / (1.0 + theta), 1.0);
    }
    std::cout << name << ": " << sampled << " sampled points; K vs orbital form " << orbital_k.ratio
              << ", DORI vs finite differences " << finite_difference.ratio << '\n';
    check(sampled >= 20, name + ": enough sampled points (got " + std::to_string(sampled) + ")");
    // The orbital form uses the MOs, the Evaluator possibly the FCHK density
    // matrices, which agree to about 1e-8 relative.
    check(orbital_k.ratio <= 1e-7,
          name + ": K agrees with -(1/2) sum phi lap phi to 1e-7 of |tau| + |lap rho| / 4 (worst " +
              num(orbital_k.ratio) + ")");
    check(finite_difference.ratio <= 1e-5,
          name + ": DORI agrees with finite differences to 1e-5 (worst " +
              num(finite_difference.ratio) + ")");

    // Identical results with 2 and 4 threads.
    for (int threads : {2, 4}) {
        kumo_test::set_threads(threads);
        const Results again = Evaluator(molecule).evaluate(grid, request);
        bool same = true;
        for (const auto& key : request.properties())
            for (std::size_t c = 0; c < r.get(key).components.size(); ++c)
                same = same &&
                       std::ranges::equal(r.get(key).component(c), again.get(key).component(c));
        check(same, name + ": identical results with " + std::to_string(threads) + " threads");
    }
    kumo_test::set_threads(1);

    // Far from the molecule (no density): the documented values.
    Grid far = grid;
    far.origin = grid.origin + Eigen::Vector3d(500.0, 0.0, 0.0);
    far.counts = Eigen::Vector3i(2, 2, 2);
    const Results empty = Evaluator(molecule).evaluate(far, request);
    bool defaults = true;
    for (std::size_t p = 0; p < 8; ++p)
        defaults = defaults && at(empty, density::total.lol(), p) == 0.0 &&
                   at(empty, density::alpha.lol(), p) == 0.0 &&
                   at(empty, density::total.interaction_region_indicator(), p) ==
                       interaction_region_indicator_cap &&
                   at(empty, density::total.density_overlap_regions_indicator(), p) == 0.0 &&
                   at(empty, density::difference.spin_polarization(), p) == 0.0 &&
                   at(empty, density::total.shannon_entropy_density(), p) == 0.0 &&
                   at(empty, density::total.fisher_information_density(), p) == 0.0 &&
                   at(empty, density::total.hamiltonian_kinetic_energy(), p) == 0.0;
    check(defaults, name + ": without density LOL, DORI, zeta, Shannon, Fisher and K are 0 and "
                           "IRI is the cap");

    // Components the fields are not defined for.
    const Evaluator evaluator(molecule);
    for (const Property& key :
         {density::difference.lol(), density::alpha.interaction_region_indicator(),
          density::beta.density_overlap_regions_indicator(),
          density::alpha.potential_energy_density(), density::difference.energy_density(),
          density::total.spin_polarization(), density::alpha.spin_polarization(),
          density::beta.shannon_entropy_density(),
          density::difference.fisher_information_density()}) {
        Request bad;
        bad.add(key);
        kumo_test::check_throws<std::invalid_argument>(
            [&] { (void)evaluator.evaluate(grid, bad); },
            name + ": a field requested for a component it is not defined for is rejected");
    }
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        for (int i = 1; i < argc; ++i)
            check_molecule(argv[i]);
        if (kumo_test::failures == 0) std::cout << "density fields OK\n";
    });
}
