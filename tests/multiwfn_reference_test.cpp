// Derived density properties and the electrostatic potential against
// Multiwfn reference Cube files (provenance in tests/REFERENCE_CUBES.md).
// The Evaluator is run on the grid read from each Cube header and compared
// point by point.
//
// The Cube files store six significant digits, so rounding alone allows a
// relative difference of 5e-6. The tolerance, 2e-5 relative, also covers the
// largest differences observed: 1.2e-5 for the reduced gradient near the O-H
// bond critical point, where |grad rho| nearly cancels, and 9e-6 for
// sign(lambda2) rho at the oxygen nucleus.
//
// The electrostatic potential also has an absolute tolerance of 1e-7
// hartree/e, for small values: the grid points are read from the Cube header,
// whose six decimals move them by up to 5e-7 bohr, and the potential changes
// by up to about 5e-8 over that distance away from the nuclei. Points within
// 1e-3 bohr of a nucleus are skipped: there the two programs treat the
// singular nuclear term differently (Multiwfn writes about 2e16 on the O2
// nucleus that its grid contains).
//
// The orbital-weighted Fukui functions are compared after removing a known
// difference of definition (see compare_orbital_weighted).
#include "chem/kumo.hpp"
#include "chem/kumo/cube.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <span>
#include <string>
#include <vector>

namespace {
using namespace chem::kumo;

constexpr double relative_tolerance = 2e-5;

struct Tolerance {
    double relative = relative_tolerance;
    double absolute = 0.0;          ///< added to relative * |expected|
    double nucleus_exclusion = 0.0; ///< bohr; points this close to a nucleus are skipped
    /// What the relative tolerance multiplies: the value, or for quantities
    /// that may be small differences of large terms, the size of those terms:
    /// |tau| + |lap rho| / 4 for K, V and H, and rho (1 + |ln rho|) for
    /// -rho ln rho, which vanishes at rho = 1, and (|ALIE| + |LEA|) / 2 for
    /// the local electronegativity (ALIE + LEA) / 2, which crosses zero.
    /// For the energy-density fields of 2.1: |G_a| + |H_aa| / 4 for the term
    /// K_a = G_a - H_aa / 4 (H_aa the diagonal of the density Hessian),
    /// (|tau| + |lap rho| / 4) / rho for H / rho, and
    /// (2 |tau| + |lap rho| / 4) / tau for |V| / G.
    enum class Scale {
        Value,
        Energy,
        Shannon,
        OrbitalEnergies,
        EnergyAxis,
        EnergyPerElectron,
        VirialRatio
    } scale = Scale::Value;
    /// With significant > 0: at least minimum_significant points must have
    /// |expected| >= significant, and there the relative tolerance alone
    /// must hold, so that an absolute term cannot hide a disagreement and the
    /// test does not compare only values near zero.
    double significant = 0.0;
    std::size_t minimum_significant = 0;
};

// `axis` selects the component of a property with several (x, y, z).
void compare(const std::filesystem::path& fchk, const std::filesystem::path& cube_path,
             Property key, const std::string& label, Tolerance tolerance = {},
             std::size_t axis = 0) {
    const auto molecule = Molecule::from_fchk(fchk);
    const auto reference = read_cube(cube_path);
    const auto& grid = reference.field.grid();
    const auto expected = reference.field.values();

    kumo_test::check(reference.field.atoms().size() == molecule.atoms().size(),
                     label + ": the Cube has the atoms of the FCHK");
    for (std::size_t a = 0; a < std::min(reference.field.atoms().size(), molecule.atoms().size());
         ++a)
        kumo_test::check(
            (reference.field.atoms()[a].position - molecule.atoms()[a].position).norm() < 1e-5,
            label + ": atom " + std::to_string(a + 1) + " is at the FCHK position");

    Request request;
    request.add(key)
        .add(density::total.value())
        .add(density::total.tau())
        .add(density::total.laplacian())
        .add(density::total.tau_components())
        .add(density::total.hessian());
    const bool orbital_scale = tolerance.scale == Tolerance::Scale::OrbitalEnergies;
    if (orbital_scale)
        request.add(density::total.average_local_ionization_energy())
            .add(density::total.local_electron_affinity());
    // Keep the Results alive: component() is a view into it.
    const Results results = Evaluator(molecule).evaluate(grid, request);
    const auto actual = results.get(key).component(axis);
    const auto rho = results.get(density::total.value()).component();
    const auto tau = results.get(density::total.tau()).component();
    const auto laplacian = results.get(density::total.laplacian()).component();
    const auto tau_axis = results.get(density::total.tau_components()).component(axis);
    // Hessian components in the order xx, xy, xz, yy, yz, zz: the diagonal
    // entry of `axis` is component 0, 3 or 5.
    const std::size_t diagonal[] = {0, 3, 5};
    const auto hessian_axis =
        results.get(density::total.hessian()).component(diagonal[std::min<std::size_t>(axis, 2)]);
    const auto alie =
        orbital_scale ? results.get(density::total.average_local_ionization_energy()).component()
                      : std::span<const double>();
    const auto lea = orbital_scale
                         ? results.get(density::total.local_electron_affinity()).component()
                         : std::span<const double>();
    if (!kumo_test::check(actual.size() == expected.size(),
                          label + ": the Evaluator returns one value per Cube point"))
        return;

    // Largest difference as a fraction of the allowed one.
    double worst = 0.0, worst_relative = 0.0, worst_significant = 0.0;
    std::size_t worst_point = 0, skipped = 0, significant = 0;
    bool same_sign = true;
    const auto& atoms = molecule.atoms();
    for (std::size_t p = 0; p < expected.size(); ++p) {
        const std::size_t k = p % static_cast<std::size_t>(grid.counts[2]);
        const std::size_t j = (p / static_cast<std::size_t>(grid.counts[2])) %
                              static_cast<std::size_t>(grid.counts[1]);
        const std::size_t i =
            p / static_cast<std::size_t>(grid.counts[2]) / static_cast<std::size_t>(grid.counts[1]);
        const Eigen::Vector3d point =
            grid.origin + grid.steps * Eigen::Vector3d(double(i), double(j), double(k));
        if (std::any_of(atoms.begin(), atoms.end(), [&](const Atom& atom) {
                return (atom.position - point).norm() < tolerance.nucleus_exclusion;
            })) {
            ++skipped;
            continue;
        }
        const double difference = std::abs(actual[p] - expected[p]);
        const double scale =
            tolerance.scale == Tolerance::Scale::Energy
                ? std::abs(tau[p]) + 0.25 * std::abs(laplacian[p])
            : tolerance.scale == Tolerance::Scale::Shannon
                ? (rho[p] > 0.0 ? rho[p] * (1.0 + std::abs(std::log(rho[p]))) : 0.0)
            : tolerance.scale == Tolerance::Scale::EnergyAxis
                ? std::abs(tau_axis[p]) + 0.25 * std::abs(hessian_axis[p])
            : tolerance.scale == Tolerance::Scale::EnergyPerElectron
                ? (rho[p] > 0.0 ? (std::abs(tau[p]) + 0.25 * std::abs(laplacian[p])) / rho[p] : 0.0)
            : tolerance.scale == Tolerance::Scale::VirialRatio
                ? (tau[p] > 0.0 ? (2.0 * std::abs(tau[p]) + 0.25 * std::abs(laplacian[p])) / tau[p]
                                : 0.0)
            : orbital_scale ? 0.5 * (std::abs(alie[p]) + std::abs(lea[p]))
                            : std::abs(expected[p]);
        const double allowed = tolerance.relative * scale + tolerance.absolute;
        const double ratio = difference / std::max(allowed, 1e-300);
        if (ratio > worst) {
            worst = ratio;
            worst_point = p;
        }
        if (std::abs(expected[p]) > 1e3 * tolerance.absolute)
            worst_relative = std::max(worst_relative, difference / std::abs(expected[p]));
        if (tolerance.significant > 0.0 && std::abs(expected[p]) >= tolerance.significant) {
            ++significant;
            worst_significant = std::max(worst_significant, difference / std::abs(expected[p]));
        }
        if (std::abs(expected[p]) > tolerance.absolute)
            same_sign = same_sign && (actual[p] < 0.0) == (expected[p] < 0.0);
    }
    std::cout << label << ": " << expected.size() - skipped << " points (" << skipped
              << " skipped near nuclei), max difference " << worst << " of the tolerance (point "
              << worst_point << ": " << actual[worst_point] << " versus " << expected[worst_point]
              << "), max relative difference " << worst_relative << '\n';
    kumo_test::check(same_sign, label + ": every value has the sign of the reference");
    if (tolerance.significant > 0.0) {
        std::cout << label << ": " << significant
                  << " points with |value| >= " << tolerance.significant
                  << ", max relative difference " << worst_significant << '\n';
        kumo_test::check(significant >= tolerance.minimum_significant,
                         label + ": at least " + std::to_string(tolerance.minimum_significant) +
                             " points with |value| >= " + kumo_test::num(tolerance.significant) +
                             " (got " + std::to_string(significant) + ")");
        kumo_test::check(worst_significant <= tolerance.relative,
                         label + ": where |value| >= " + kumo_test::num(tolerance.significant) +
                             ", agrees with Multiwfn within " + kumo_test::num(tolerance.relative) +
                             " relative alone (max " + kumo_test::num(worst_significant) + ")");
    }
    kumo_test::check(worst <= 1.0, label + ": agrees with Multiwfn within " +
                                       kumo_test::num(tolerance.relative) + " relative + " +
                                       kumo_test::num(tolerance.absolute) + " absolute (max " +
                                       kumo_test::num(worst) + " of the tolerance)");
}

// Orbital-weighted Fukui functions f+ and f- (user functions 95 and 96)
// and, if `all`, f0 and the dual descriptor (97 and 98), with Multiwfn's
// default delta of 0.1 hartree.
//
// kumo sums over every orbital; Multiwfn keeps only the orbitals whose
// unnormalized Gaussian t_i = exp(-((mu - e_i) / delta)^2) is at least 1e-3,
// scanning up from the LUMO (f+) or down from the HOMO (f-) and stopping at
// the first one below, and normalizes the weights over the kept orbitals.
// The test rebuilds that truncated sum from the orbitals evaluated by kumo,
//   m+ = sum_kept t_i |phi_i|^2 / sum_kept t_i   (likewise m-),
// and adds the difference d+ = f+ - m+ (and d-) to the Cube values: then
// the Cube value plus d+ must equal kumo's f+ within the Cube rounding.
// Without this step the two programs differ by the truncated weight, 0.4 %
// of f+ for acrolein, far above the Cube precision.
//
// Where |value| >= 1e-4 bohr^-3 (more than 6000 points per Cube, the largest
// values being 0.14 to 0.47) the relative tolerance 2e-5 must hold alone.
// Below, at a few points with values of 1e-6 to 7e-5, the programs differ by
// up to 3e-9 bohr^-3 absolute, 4e-5 relative; an absolute term of 1e-8
// bohr^-3, about three times that and 2e-8 to 7e-8 of the largest values,
// covers them.
void compare_orbital_weighted(const std::filesystem::path& fchk,
                              const std::filesystem::path& directory, const std::string& prefix,
                              const std::string& name, bool all) {
    const auto molecule = Molecule::from_fchk(fchk);
    const auto& energies = molecule.orbital_energies(SpinChannel::Alpha);
    const std::size_t homo = molecule.occupied_count(SpinChannel::Alpha); // one-based
    const std::size_t count = molecule.orbital_count(SpinChannel::Alpha);
    constexpr double delta = 0.1, cut = 1e-3;
    constexpr double absolute = 1e-8, significant_value = 1e-4;
    constexpr std::size_t minimum_significant = 5000;
    const double mu = 0.5 * (energies[homo - 1] + energies[homo]);
    const auto gaussian = [&](std::size_t index) { // one-based
        const double x = (mu - energies[index - 1]) / delta;
        return std::exp(-x * x);
    };
    std::vector<std::size_t> plus_kept, minus_kept;
    for (std::size_t i = homo + 1; i <= count && gaussian(i) >= cut; ++i)
        plus_kept.push_back(i);
    for (std::size_t i = homo; i >= 1 && gaussian(i) >= cut; --i)
        minus_kept.push_back(i);
    if (!kumo_test::check(!plus_kept.empty() && !minus_kept.empty(),
                          name + ": Multiwfn keeps at least the LUMO and the HOMO"))
        return;
    std::cout << name << ": Multiwfn keeps " << plus_kept.size() << " orbitals in f+ and "
              << minus_kept.size() << " in f-\n";

    struct Function {
        const char* file;
        Property key;
        const char* label;
    };

    const auto& t = density::total;
    std::vector<Function> functions = {{"ow_fukui_plus", t.orbital_weighted_fukui_plus(), "f+"},
                                       {"ow_fukui_minus", t.orbital_weighted_fukui_minus(), "f-"}};
    if (all) {
        functions.push_back({"ow_fukui_zero", t.orbital_weighted_fukui_zero(), "f0"});
        functions.push_back(
            {"ow_dual_descriptor", t.orbital_weighted_dual_descriptor(), "dual descriptor"});
    }
    for (const auto& function : functions) {
        const auto label = name + " orbital-weighted " + function.label;
        const auto reference = read_cube(directory / (prefix + "_" + function.file + ".cube.gz"));
        const auto& grid = reference.field.grid();
        const auto expected = reference.field.values();
        Request request;
        request.add(function.key)
            .add(t.orbital_weighted_fukui_plus())
            .add(t.orbital_weighted_fukui_minus());
        for (const auto i : plus_kept)
            request.add(orbital::alpha(i).value());
        for (const auto i : minus_kept)
            request.add(orbital::alpha(i).value());
        // Keep the Results alive: component() is a view into it.
        const Results results = Evaluator(molecule).evaluate(grid, request);
        const auto actual = results.get(function.key).component();
        const auto plus = results.get(t.orbital_weighted_fukui_plus()).component();
        const auto minus = results.get(t.orbital_weighted_fukui_minus()).component();
        if (!kumo_test::check(actual.size() == expected.size(),
                              label + ": the Evaluator returns one value per Cube point"))
            return;
        const auto truncated = [&](const std::vector<std::size_t>& kept, std::size_t p) {
            double sum = 0.0, norm = 0.0;
            for (const auto i : kept) {
                const double phi = results.get(orbital::alpha(i).value()).component()[p];
                sum += gaussian(i) * phi * phi;
                norm += gaussian(i);
            }
            return sum / norm;
        };
        double worst = 0.0, worst_significant = 0.0, worst_correction = 0.0;
        std::size_t significant = 0;
        for (std::size_t p = 0; p < expected.size(); ++p) {
            const double d_plus = plus[p] - truncated(plus_kept, p);
            const double d_minus = minus[p] - truncated(minus_kept, p);
            const double correction = function.key == t.orbital_weighted_fukui_plus()    ? d_plus
                                      : function.key == t.orbital_weighted_fukui_minus() ? d_minus
                                      : function.key == t.orbital_weighted_fukui_zero()
                                          ? 0.5 * (d_plus + d_minus)
                                          : d_plus - d_minus;
            // f+, f- and f0 are positive and compared relative to their
            // value; the dual descriptor crosses zero and is compared on the
            // scale |f+| + |f-|.
            const double scale = function.key == t.orbital_weighted_dual_descriptor()
                                     ? std::abs(plus[p]) + std::abs(minus[p])
                                     : std::abs(actual[p]);
            const double difference = std::abs(actual[p] - (expected[p] + correction));
            const double relative = difference / std::max(relative_tolerance * scale, 1e-300);
            worst = std::max(worst, difference / (relative_tolerance * scale + absolute));
            if (std::abs(expected[p]) >= significant_value) {
                ++significant;
                worst_significant = std::max(worst_significant, relative);
            }
            worst_correction =
                std::max(worst_correction, std::abs(correction) / std::max(scale, 1e-300));
        }
        std::cout << label << ": " << expected.size() << " points, max difference " << worst
                  << " of the tolerance; " << significant
                  << " points with |value| >= " << significant_value << ", max difference there "
                  << worst_significant
                  << " of the relative tolerance alone; largest truncation correction "
                  << worst_correction << " of the value\n";
        kumo_test::check(significant >= minimum_significant,
                         label + ": at least " + std::to_string(minimum_significant) +
                             " points with |value| >= " + kumo_test::num(significant_value) +
                             " (got " + std::to_string(significant) + ")");
        kumo_test::check(worst_significant <= 1.0,
                         label + ": where |value| >= " + kumo_test::num(significant_value) +
                             ", agrees with Multiwfn within " + kumo_test::num(relative_tolerance) +
                             " relative alone after the truncation correction (max " +
                             kumo_test::num(worst_significant) + " of the tolerance)");
        kumo_test::check(worst <= 1.0, label + ": agrees with Multiwfn within " +
                                           kumo_test::num(relative_tolerance) + " relative + " +
                                           kumo_test::num(absolute) +
                                           " absolute after the truncation correction (max " +
                                           kumo_test::num(worst) + " of the tolerance)");
    }
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 4, 4,
                                 "WATER.fchk O2.fchk ACROLEIN.fchk MULTIWFN_DIRECTORY"))
        return 2;
    return kumo_test::run([&] {
        const std::filesystem::path water = argv[1], o2 = argv[2], acrolein = argv[3],
                                    directory = argv[4];
        compare(water, directory / "h2o_elf.cube.gz", density::total.elf(), "water ELF");
        compare(water, directory / "h2o_rdg.cube.gz", density::total.reduced_gradient(),
                "water reduced gradient");
        compare(water, directory / "h2o_sign_lambda2_rho.cube.gz",
                density::total.sign_lambda2_density(), "water sign(lambda2) rho");
        // Open shell: the total ELF uses the spin-polarized form.
        compare(o2, directory / "o2_elf.cube.gz", density::total.elf(), "O2 ELF");
        const Tolerance potential{relative_tolerance, 1e-7, 1e-3};
        compare(water, directory / "h2o_esp.cube.gz", density::total.electrostatic_potential(),
                "water electrostatic potential", potential);
        compare(o2, directory / "o2_esp.cube.gz", density::total.electrostatic_potential(),
                "O2 electrostatic potential", potential);
        // DORI is in [0, 1]; theta combines the Hessian and the gradient with
        // cancellation, so an absolute term covers small values.
        const Tolerance dori{relative_tolerance, 1e-5};
        // K, V and H combine tau and lap rho / 4 with cancellation.
        Tolerance energy{5e-5};
        energy.scale = Tolerance::Scale::Energy;
        Tolerance shannon;
        shannon.scale = Tolerance::Scale::Shannon;
        // As the reduced gradient, largest near the O-H bond critical point.
        const Tolerance iri{5e-5};
        // Spin polarization in [-1, 1]: a small ratio of nearly equal alpha
        // and beta densities, compared with an absolute term.
        const Tolerance polarization{relative_tolerance, 1e-6};
        // |grad rho|^2 doubles the relative error of the gradient near the
        // O-H bond critical point (1.2e-5 for RDG).
        const Tolerance fisher{1e-4};
        compare(water, directory / "h2o_lol.cube.gz", density::total.lol(), "water LOL");
        compare(o2, directory / "o2_lol.cube.gz", density::total.lol(), "O2 LOL");
        compare(water, directory / "h2o_iri.cube.gz", density::total.interaction_region_indicator(),
                "water IRI", iri);
        compare(water, directory / "h2o_dori.cube.gz",
                density::total.density_overlap_regions_indicator(), "water DORI", dori);
        compare(water, directory / "h2o_kinetic_k.cube.gz",
                density::total.hamiltonian_kinetic_energy(), "water K", energy);
        compare(o2, directory / "o2_kinetic_k.cube.gz", density::total.hamiltonian_kinetic_energy(),
                "O2 K", energy);
        compare(water, directory / "h2o_potential_energy.cube.gz",
                density::total.potential_energy_density(), "water V", energy);
        compare(water, directory / "h2o_energy_density.cube.gz", density::total.energy_density(),
                "water H", energy);
        compare(o2, directory / "o2_spin_polarization.cube.gz",
                density::difference.spin_polarization(), "O2 spin polarization", polarization);
        compare(water, directory / "h2o_shannon.cube.gz", density::total.shannon_entropy_density(),
                "water Shannon entropy density", shannon);
        compare(water, directory / "h2o_fisher.cube.gz",
                density::total.fisher_information_density(), "water Fisher information density",
                fisher);
        // Orbital-energy fields (restricted acrolein, unrestricted O2). ALIE,
        // LEA and the local hardness, of order 1 hartree, take the default
        // tolerance. LEAE vanishes on the nodal planes of the negative-energy
        // virtual orbitals, where both programs give rounding noise (-4e-39
        // versus -1e-34, at most 2.5e-33 apart), so it adds 1e-10 hartree
        // absolute, 1e-8 of its largest magnitude (0.014 hartree in
        // acrolein, 0.024 in O2). Most points carry a substantial LEAE: where
        // |LEAE| >= 1e-3 hartree (12258 acrolein and 11788 O2 points) the
        // relative tolerance must hold alone, and the largest difference
        // there is 6.4e-6 relative, the Cube rounding. The local
        // electronegativity (ALIE + LEA) / 2 crosses zero and is compared on
        // the scale of its terms, (|ALIE| + |LEA|) / 2 (largest difference
        // 2.7e-6 of it).
        const Tolerance orbital{};
        Tolerance attachment{relative_tolerance, 1e-10};
        attachment.significant = 1e-3;
        attachment.minimum_significant = 10000;
        Tolerance electronegativity;
        electronegativity.scale = Tolerance::Scale::OrbitalEnergies;
        compare(acrolein, directory / "acrolein_alie.cube.gz",
                density::total.average_local_ionization_energy(), "acrolein ALIE", orbital);
        compare(o2, directory / "o2_alie.cube.gz", density::total.average_local_ionization_energy(),
                "O2 ALIE", orbital);
        compare(acrolein, directory / "acrolein_lea.cube.gz",
                density::total.local_electron_affinity(), "acrolein LEA", orbital);
        compare(o2, directory / "o2_lea.cube.gz", density::total.local_electron_affinity(),
                "O2 LEA", orbital);
        compare(acrolein, directory / "acrolein_leae.cube.gz",
                density::total.local_electron_attachment_energy(), "acrolein LEAE", attachment);
        compare(o2, directory / "o2_leae.cube.gz",
                density::total.local_electron_attachment_energy(), "O2 LEAE", attachment);
        compare(acrolein, directory / "acrolein_local_electronegativity.cube.gz",
                density::total.local_electronegativity(), "acrolein local electronegativity",
                electronegativity);
        compare(acrolein, directory / "acrolein_local_hardness.cube.gz",
                density::total.local_hardness(), "acrolein local hardness", orbital);
        // Tsirelson-Stash ELF and LOL (closed shells; Multiwfn's ELFLOL_type= 1).
        compare(water, directory / "h2o_elf_ts.cube.gz", density::total.elf_tsirelson_stash(),
                "water Tsirelson-Stash ELF");
        compare(acrolein, directory / "acrolein_elf_ts.cube.gz",
                density::total.elf_tsirelson_stash(), "acrolein Tsirelson-Stash ELF");
        compare(water, directory / "h2o_lol_ts.cube.gz", density::total.lol_tsirelson_stash(),
                "water Tsirelson-Stash LOL");
        compare(acrolein, directory / "acrolein_lol_ts.cube.gz",
                density::total.lol_tsirelson_stash(), "acrolein Tsirelson-Stash LOL");
        // Energy densities of 2.1 (Multiwfn function 7 and user functions
        // 81-86, 35, 17, 13).
        // The Cartesian terms of G reach 1.7e-5 relative (water G_x, where
        // it is 0.02 hartree/bohr^3): 5e-5, as for K, V and H. K_a, H / rho
        // and |V| / G cancel terms as K, V and H do and are compared on the
        // scale of those terms, with the same 5e-5.
        const Tolerance energy_term{5e-5};
        Tolerance energy_axis{5e-5};
        energy_axis.scale = Tolerance::Scale::EnergyAxis;
        Tolerance per_electron{5e-5};
        per_electron.scale = Tolerance::Scale::EnergyPerElectron;
        Tolerance virial{5e-5};
        virial.scale = Tolerance::Scale::VirialRatio;
        compare(water, directory / "h2o_kinetic_g.cube.gz", density::total.tau(), "water G");
        // The O2 grid contains one nucleus. The six decimals of the Cube
        // header move it by up to 5e-7 bohr, and on the nuclear cusp G
        // changes by about 2Z = 32 relative per bohr: 2.9e-5 at the nucleus,
        // above the tolerance. That point is skipped (as for the ESP); every
        // other point agrees within 2e-5.
        compare(o2, directory / "o2_kinetic_g.cube.gz", density::total.tau(), "O2 G",
                Tolerance{relative_tolerance, 0.0, 1e-3});
        const char* axes[] = {"x", "y", "z"};
        for (std::size_t a = 0; a < 3; ++a) {
            compare(water, directory / (std::string("h2o_kinetic_g_") + axes[a] + ".cube.gz"),
                    density::total.tau_components(), std::string("water G_") + axes[a], energy_term,
                    a);
            compare(water, directory / (std::string("h2o_kinetic_k_") + axes[a] + ".cube.gz"),
                    density::total.hamiltonian_kinetic_energy_components(),
                    std::string("water K_") + axes[a], energy_axis, a);
        }
        compare(o2, directory / "o2_kinetic_k_x.cube.gz",
                density::total.hamiltonian_kinetic_energy_components(), "O2 K_x", energy_axis, 0);
        compare(water, directory / "h2o_potential_kinetic_ratio.cube.gz",
                density::total.potential_kinetic_ratio(), "water |V| / G", virial);
        compare(o2, directory / "o2_potential_kinetic_ratio.cube.gz",
                density::total.potential_kinetic_ratio(), "O2 |V| / G", virial);
        compare(water, directory / "h2o_energy_density_per_electron.cube.gz",
                density::total.energy_density_per_electron(), "water H / rho", per_electron);
        compare(water, directory / "h2o_kinetic_g_per_electron.cube.gz",
                density::total.kinetic_energy_per_electron(), "water G / rho");
        // Orbital-weighted Fukui functions (closed shell only).
        compare_orbital_weighted(acrolein, directory, "acrolein", "acrolein", true);
        compare_orbital_weighted(water, directory, "h2o", "water", false);
    });
}
