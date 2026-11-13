# chem-kumo C++ API guide

Include `<chem/kumo.hpp>` and link against `kumo_dep` in Meson. Public types
are in namespace `chem::kumo`. Examples use fully qualified names (such as
`chem::kumo::Molecule` and `chem::kumo::Spin`) rather than namespace-wide
`using` directives. All spatial coordinates and grid vectors use atomic
units (bohr).

## Reading an FCHK

`Molecule::from_fchk(path)` accepts plain FCHK or gzip-compressed FCHK (detected by magic bytes), and constructs the molecule, basis evaluator, SCF
matrices and available MO coefficient/energy arrays. `FormattedCheckpoint`
provides typed access to named scalar and array records; a missing or
incorrectly typed record raises an exception. No Gaussian executable is
required at runtime.

## Grid and fields

`Grid::around(atoms, margin_bohr, spacing_bohr)` builds an axis-aligned
regular grid around the atomic bounding box. `Grid::origin` is the first
point; `Grid::steps.col(k)` is the increment along grid axis `k`;
`Grid::counts` gives the point count along each axis.

A grid need not be axis-aligned or three-dimensional: an axis with one point
gives a plane, two give a line. Two constructors build them (1.13.0):

- `Grid::line(start, end, margin_bohr, spacing_bohr)`: points on the line
  through `start` and `end`, from `margin_bohr` before `start` to
  `margin_bohr` after `end`, evenly spaced with both ends included; the
  fewest points (at least 2) whose spacing does not exceed `spacing_bohr`.
  Counts are (n, 1, 1); the first axis runs from start to end and the other
  two are orthogonal unit vectors.
- `Grid::plane(a, b, c, margin_bohr, spacing_bohr)`: a rectangle in the
  plane of three points. The first axis runs along b - a, the second
  perpendicular to it towards c, and the third is the unit normal u x v,
  with one point. The rectangle covers the three points enlarged by
  `margin_bohr` on every side, edges included, with the fewest points along
  each axis whose spacing does not exceed `spacing_bohr`. Counts are
  (nu, nv, 1), so values are stored with v fastest.

Both throw `std::invalid_argument` for nonfinite points, coincident ends or
(nearly) collinear points, a negative or nonfinite margin, a nonpositive
spacing, or more points than an `int` count holds. A length that is a
multiple of the spacing up to rounding keeps that multiple, so a line or
plane through points of another regular grid reproduces them. Every
`Evaluator` property can be evaluated on these grids (for a single point,
use a `Grid` with `origin` at the point and the default steps and counts);
`read_cube()` and `FieldSampler` need at least two points per axis. The test
`line-plane` checks the geometry, agreement with point-by-point evaluation
(1e-12, and 2e-8 for the electrostatic potential), and LOL and ELF on a
line and a plane through the points of the Multiwfn water Cube files
(largest difference 24% of the 2e-5 relative tolerance).

`Field::operator()(i,j,k)` accesses a sampled value. `Field::values()`
provides a read-only contiguous view. `write_cube(path)` writes scalar fields;
`write_orbital_cube(path, orbital_index)` writes a Gaussian-style orbital
Cube with an MO index record. `Field` stores values, atoms, and grid geometry.

## Density methods

`Spin::Total`, `Spin::Alpha`, `Spin::Beta`, and `Spin::Difference` select
`rho_alpha + rho_beta`, `rho_alpha`, `rho_beta`, and
`rho_alpha - rho_beta`, respectively. Densities have units e/bohr^3.

- `density(grid, spin)` returns one component.
- `densities(grid, {components...})` computes selected fields with shared AO
  evaluation; duplicate component requests do not create duplicate fields.
- `spin_densities(grid)` returns all four components in a named structure.
- `densities_blocked(grid, components, block_size, screening_threshold,
  progress)` explicitly selects a block size (the default density path uses
  128 points per block), an AO screening threshold and an optional progress
  callback.
- Progress reports: `ProgressCallback` (`<chem/kumo/density.hpp>`) is a
  `std::function<void(std::size_t done, std::size_t total)>`, accepted by
  `densities_blocked`, `orbitals` and `EvaluationOptions::progress`. It is
  called after each block of grid points with the number of points done so
  far and the total; calls are serialized (never concurrent, possibly from
  an OpenMP worker thread), `done` increases up to `total`, and an
  exception thrown by the callback aborts the evaluation and reaches the
  caller. For `Evaluator::evaluate`, `total` counts the points of all passes
  over the grid: the derivative or AO-MO pass, if any, plus one pass per
  AO-P density or orbital field. The callback does not change the
  results.
- AO screening: in each block, AOs whose absolute value is at most
  `default_screening_threshold` (1e-12) at every point are left out of the
  density-matrix contraction. This makes large-basis densities several times
  faster; the change in the density is at rounding level (at most about
  1e-13 e/bohr^3 in the tested systems). A threshold of 0 gives the exact
  contraction. `EvaluationOptions::screening_threshold` sets the same
  threshold for `Evaluator` AO_P densities and density derivatives, and
  `density_derivatives_blocked()` takes it as a parameter. For derivatives
  an AO is left out only if its value and all its derivatives are at most
  the threshold in the block; on HBC def2-TZVP this makes AO_P gradients and
  Hessians about five times faster, with relative changes of about 1e-15.
  The AO_MO contraction (densities and derivatives) is screened the same
  way: negligible AOs are dropped from the AO panels and from the occupied
  MO coefficients of each block (about 1.6 times faster on HBC
  def2-TZVP, where AO evaluation then dominates).
- `density_pointwise(grid, spin)` retains the original pointwise algorithm
  for regression comparisons.

## Molecular orbitals

`orbital_count(channel)`, `occupied_count(channel)`, and
`orbital_energies(channel)` inspect alpha or beta MO metadata.
`orbital(grid, index_one_based, channel)` returns a signed spatial orbital
amplitude, and `orbitals(grid, indices_one_based, channel)` computes multiple
orbital fields with shared AO batches. Index 1 means the first orbital.
`Spin::Alpha` and `Spin::Beta` are valid MO channels; `Spin::Total` and
`Spin::Difference` describe densities and are not MO channels.

An orbital is `psi_i(r) = sum_mu C[mu,i] chi_mu(r)`. The returned value is
**not** its square and has units bohr^(-3/2). No phase adjustment or
normalization is applied. For restricted calculations alpha and beta refer
to the same spatial orbital set; unrestricted calculations use independent
alpha/beta coefficients and orbital energies. HOMO/LUMO labels in the CLI
are determined from occupied counts per channel.

## AO and basis interfaces

`BasisSet` represents Gaussian shells, including Cartesian, pure, and SP
shell metadata. `BasisEvaluator::evaluate(point)` returns AO values in FCHK
ordering; `evaluate_into` fills a caller-provided vector; `evaluate_batch`
produces AO columns for a batch of points. The `evaluate_into_reference`
method is a numerical-regression implementation.

Shell cutoff (1.6.0): for each shell, the evaluator precomputes a radius
beyond which a rigorous bound on every AO of the shell and on its first and
second derivatives is at most `BasisEvaluator::ao_cutoff` (1e-14). At points
farther than that from the shell center, `evaluate*` and
`evaluate_derivatives_batch` return exact zeros for the shell instead of
computing them; `evaluate_into_reference` does not apply the cutoff. On HBC
def2-TZVP grids this makes densities about 1.3 times (AO_MO about 1.5
times) and derivatives 1.1 to 1.3 times faster, with changes of at most
about 6e-14 in the results.

## Spatial derivative API (1.2.0)

`BasisEvaluator::evaluate_derivatives(point)` returns an `AODerivatives`
structure: `values` (N AO), `gradient` (N AO x 3) and `hessian`
(N AO x 6). The Hessian columns are **xx, yy, zz, xy, xz, yz**; the
coordinate axes and units are Cartesian bohr. AO ordering follows the FCHK.

`BasisEvaluator::evaluate_derivatives_batch(points, order, output)` computes
the same quantities for a 3 x N batch of points without per-point
allocations, which is how the library evaluates density and orbital
derivatives on grids. `order` is 1 (values and gradient) or 2 (also the
Hessian). `output` must be N AO x (panels * N), with
`panels = BasisEvaluator::derivative_panel_count(order)` (4 or 10); it is
filled with consecutive N AO x N panels in the order value, x, y, z, xx, yy,
zz, xy, xz, yz, so component `c` at point `j` is column `c * N + j`.

`Molecule::orbital_derivatives(point, index_one_based, channel)` returns
`OrbitalDerivatives` with signed `value`, `gradient` (3-vector) and
`hessian` (symmetric 3 x 3 matrix). The MO channel must be Alpha or Beta.

`Molecule::density_derivatives(point, component)` returns
`DensityDerivatives` with `value`, `gradient`, `hessian`, `laplacian`,
and `tau`. `laplacian` equals the trace of `hessian`. The positive
kinetic-energy density is defined as
`tau = (1/2) sum_i occupation_i |grad psi_i|^2` for a physical spin
channel; for `Spin::Difference`, the returned tau is the **signed**
difference `tau_alpha - tau_beta`. Units are density bohr^-3, gradient
bohr^-4, Hessian/Laplacian bohr^-5 and tau hartree/bohr^3 in atomic units.

`Molecule::density_derivatives_blocked(points, component, block_size=128,
screening_threshold=default_screening_threshold)`
accepts a **3 x N** Eigen matrix with one Cartesian point per column,
and returns `std::vector<DensityDerivatives>` in the same order. A zero
block size, non-finite points and a negative or nonfinite screening
threshold are rejected. The method currently
materializes the result vector; it is not a streaming grid writer.

`Molecule::orbital_derivatives_blocked(points, index_one_based, channel,
block_size=128)` accepts the same 3 x N point layout and returns ordered
`OrbitalDerivatives` for one alpha/beta MO, evaluating the AO derivatives of
each block with `evaluate_derivatives_batch()`.

When built with `-Dopenmp=true`, `densities_blocked()`,
`density_derivatives_blocked()`, `orbitals()` (and `orbital()` through it), and
`orbital_derivatives_blocked()` schedule independent grid/point blocks across
OpenMP workers. OpenMP is optional and off by default. Pointwise evaluators
remain serial; thread count is controlled through `OMP_NUM_THREADS`.
See [docs/OPENMP.md](OPENMP.md) for execution and deployment notes.

The blocked derivative path uses Eigen to apply the density matrix once
to a packed panel containing AO values and their three first derivatives.
Six AO second-derivative panels are used in linear reductions; no density
matrix product is needed for those panels. AO derivatives of a block are
evaluated together with `evaluate_derivatives_batch()`, and negligible AOs
are screened out of the contraction as described above. There is no sparse
matrix contraction or new numerical dependency.

These are **spatial derivatives at fixed geometry and fixed wavefunction**,
not derivatives with respect to nuclear positions. The `kumo` CLI writes
the gradient norm, Laplacian and `tau` of any density component as Cube
files with `--property` (see the README).
The mathematical derivation is in
[chem-kumo-mathematics.tex](mathematics/chem-kumo-mathematics.tex).
Performance methodology and HBC measurements are in
[SPATIAL_DERIVATIVES.md](../benchmarks/SPATIAL_DERIVATIVES.md).

## Evaluator (1.4.0)

`Evaluator` (`<chem/kumo/evaluator.hpp>`) evaluates typed requests on a
`Grid` and is the interface the CLI uses for all density output.

```cpp
Evaluator evaluator(molecule, {.potential_tolerance = 1e-6});
Request request;
request.add(density::total.value())
    .add(density::total.gradient())
    .add(density::total.electrostatic_potential())
    .add(orbital::alpha(5).value());
const Results results = evaluator.evaluate(grid, request);
const std::span<const double> rho = results.get(density::total.value()).component();
```

- **Requests.** `density::total`, `alpha`, `beta` and `spin` (alpha minus
  beta) give `value()`, `gradient()`, `hessian()`, `laplacian()`, `tau()`,
  `elf()`, `reduced_gradient()`, `sign_lambda2_density()`,
  `electrostatic_potential()` and the fields of the sections below; `orbital::alpha(i)` and `orbital::beta(i)`
  (one-based) give `value()`. A property requested twice is evaluated once,
  and the order of a request does not matter.
- **Results.** `Results::get(property)` returns a `PropertyField`, one array
  per component in Field order (third grid axis fastest): gradients have 3
  components (x, y, z), Hessians 6 (xx, xy, xz, yy, yz, zz), the others 1.
  `component()` returns a `std::span` into the Results, which must stay
  alive while the span is used. `Results::backend_used()` reports the
  density contraction actually used.
- **Options** (`EvaluationOptions`): `backend` (`AUTO`, `AO_P`, `AO_MO`; see
  the README), `block_size` (grid points per block, default 128),
  `screening_threshold` (AO screening, default 1e-12, 0 exact),
  `potential_tolerance` (electrostatic potential error bound, default 1e-8
  hartree/e, 0 exact), `progress` (see Density methods) and `threads`
  (only 0 and 1 are accepted; use `OMP_NUM_THREADS` in OpenMP builds).
- **Lifetime.** The Evaluator keeps a pointer to the `Molecule`, which must
  outlive it and every call to `evaluate()`.
- **Errors.** The constructor throws `std::invalid_argument` for invalid
  options or for `AO_MO` on a molecule whose orbitals are not AO-MO
  eligible (`Molecule::electronic_structure_info()` reports the check: its
  differences, relative to the norm of the total density matrix, must be
  at most `orbital_check_relative_tolerance`, 1e-8, for both the total and
  the spin density matrix); `evaluate()` throws it for an invalid grid, an orbital index of
  0, or a derived property requested for a component it is not defined for.
- **Determinism.** Results do not depend on the number of OpenMP threads,
  and the electrostatic potential does not depend on the block size either
  (see [OPENMP.md](OPENMP.md)). With AO screening, different requests may differ
  at rounding level, because a block keeps the AOs needed by the most
  demanding property requested.
- **Arbitrary points (1.13.0).** `evaluate(points, request)` evaluates the
  same requests at the columns of a `3 x N` matrix (bohr), in column order,
  for example the vertices of a surface. Values equal grid evaluation at the
  same points to rounding (the test `point-evaluation` finds at most 5e-15
  relative, with both backends): densities come from the point kernels
  (with AO-P too), orbitals from their own pass. The electrostatic potential
  groups the points by cubes of `potential_cluster_edge` (1.6 bohr) aligned
  with the axes at multiples of the edge, each evaluated around its centre
  within `potential_tolerance`; it agrees with grid evaluation within 2e-10
  in that test. A nonfinite point throws `std::invalid_argument`; an empty
  matrix gives empty fields. Results do not depend on the number of OpenMP
  threads.

## Derived properties: ELF, reduced gradient, sign(lambda2) rho (1.8.0)

The Evaluator computes three derived quantities pointwise from the density
derivatives it already evaluates, with the same backends, screening and
OpenMP scheduling (results do not depend on the number of threads). With
`tau = (1/2) sum_i |grad psi_i|^2` as above:

- `density::total.elf()`, `density::alpha.elf()`, `density::beta.elf()`
  (`PropertyKind::Elf`): the electron localization function of Becke and
  Edgecombe (J. Chem. Phys. 92, 5397 (1990)), dimensionless, in [0, 1].
  Each spin channel s with `rho_s > 0` contributes the Pauli kinetic energy
  density `D_s = tau_s - |grad rho_s|^2 / (8 rho_s)` and the uniform-gas
  value `D0_s = (3/10) (6 pi^2)^(2/3) rho_s^(5/3)`, and
  `ELF = 1 / (1 + (D / D0)^2)`. For alpha or beta, D and D0 are those of the
  channel. For the total ELF they are summed over both channels
  (spin-polarized form). For a closed shell it
  reduces to `D = tau - |grad rho|^2 / (8 rho)`,
  `D0 = (3/10) (3 pi^2)^(2/3) rho^(5/3)`. Channels without density
  contribute nothing; where there is no density at all, ELF is 0. ELF is not
  defined for the spin difference.
- `density::total.reduced_gradient()` (`PropertyKind::ReducedGradient`): the
  reduced density gradient `s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3))`
  of the total density, dimensionless (Johnson et al., J. Am. Chem. Soc. 132,
  6498 (2010)). Values above `reduced_gradient_cap` (100), and points where
  `rho <= 0`, where s is undefined, are reported as 100. NCI isovalues are
  about 0.3 to 0.6, so the cap does not affect NCI isosurfaces, and Cube
  files stay finite.
- `density::total.sign_lambda2_density()`
  (`PropertyKind::SignLambda2Density`): `sign(lambda2) * rho` of the total
  density, in e/bohr^3, where `lambda1 <= lambda2 <= lambda3` are the
  eigenvalues of the density Hessian. `lambda2 == 0` counts as positive.
  Negative values mark bonding (attractive) interactions in NCI analysis.

The reduced gradient and `sign(lambda2) rho` are defined only for the total
density. `Evaluator::evaluate` throws `std::invalid_argument` if a derived
property is requested for a component it is not defined for.

Accuracy. On the grids of the tests, the derived values equal their formulas
applied to the Evaluator's own value, gradient, Hessian and tau to about
1e-12. The Pauli term agrees, at points with `rho_s >= 1e-3`, with the
cancellation-free orbital form
`(1 / (2 rho_s)) sum_{i<j} |psi_i grad psi_j - psi_j grad psi_i|^2` to
better than 1e-9 relative to tau. With the default AO screening, ELF differs
from the unscreened value by at most 2e-8, and s by at most 5e-10 relative,
on the default `kumo` grids of O2 and acrolein. These grids include points
with densities down to 1e-18. With AO_P, densities requested alongside come
from a separate grid kernel and agree with the values used internally for
the derived properties only to rounding.

On water (B3LYP/def2-TZVP), the O-H bond critical point is found where s
nearly vanishes, with `sign(lambda2) rho` = -0.361 e/bohr^3 and ELF 0.98.

The test `multiwfn-reference` compares ELF, s and `sign(lambda2) rho` with
Multiwfn 2026.9.1 Cube files of water and O2 (open shell), point by point.
They agree within 1.2e-5 relative, at the six significant digits of the
reference files (see [tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)).

## Density fields: LOL, IRI, DORI, energy densities, spin polarization (1.12.0)

Nine more fields are computed pointwise from the density derivatives, in
the same pass as ELF and with the same backends, screening and OpenMP
scheduling (results do not depend on the number of threads). They follow
Multiwfn's definitions, which the reference tests check:

- `lol()` (`PropertyKind::Lol`), for Total, Alpha or Beta: the localized
  orbital locator of Schmider and Becke (J. Mol. Struct. THEOCHEM 527, 51
  (2000)), `LOL = t / (1 + t)` with `t = D0 / tau`, in [0, 1]. D0 is the
  uniform-gas value of ELF and, as for ELF, the total LOL sums tau and D0
  over both channels (spin-polarized form; for a closed shell the familiar
  `D0 = (3/10) (3 pi^2)^(2/3) rho^(5/3)`). It is 0 where there is no
  density and 1 where tau = 0 but the density is positive (the limit;
  Multiwfn returns 0 there). Not defined for the spin difference.
- `interaction_region_indicator()` (`PropertyKind::InteractionRegionIndicator`,
  Total only): `IRI = |grad rho| / rho^1.1` (Lu and Chen, Chemistry-Methods
  1, 231 (2021)), with isosurfaces of about 1 showing covalent bonds and
  weak interactions together. Values above `interaction_region_indicator_cap`
  (100), and points where `rho <= 0`, are reported as 100. Multiwfn sets IRI
  to 5 where rho is below `IRI_rhocut` (5e-5 by default); kumo applies no
  such cut.
- `density_overlap_regions_indicator()`
  (`PropertyKind::DensityOverlapRegionsIndicator`, Total only): DORI of de
  Silva and Corminboeuf (J. Chem. Theory Comput. 10, 3745 (2014)),
  `theta / (1 + theta)` with
  `theta = |grad (|grad rho|^2 / rho^2)|^2 / (|grad rho|^2 / rho^2)^3`, in
  [0, 1]. It is computed from `q = grad rho / rho` as
  `theta = 4 |(H / rho) q - |q|^2 q|^2 / |q|^6`, which does not depend on
  the scale of rho. At a critical point of rho (q = 0) DORI is 1, its limit;
  where `rho <= 0` it is 0, its limit far from the molecule.
- `hamiltonian_kinetic_energy()` (`PropertyKind::HamiltonianKineticEnergy`,
  any component): the kinetic energy density
  `K = -(1/2) sum_i psi_i lap psi_i = tau - lap rho / 4` (hartree/bohr^3),
  which, unlike tau (G in Bader's notation), can be negative.
- `potential_energy_density()` (`PropertyKind::PotentialEnergyDensity`,
  Total only): `V = lap rho / 4 - 2 tau`, from the local virial theorem.
- `energy_density()` (`PropertyKind::EnergyDensity`, Total only): the
  electronic energy density `H = G + V = lap rho / 4 - tau = -K`. Its sign
  at a bond critical point is used to characterize interactions (Cremer and
  Kraka).
- `density::spin.spin_polarization()` (`PropertyKind::SpinPolarization`,
  Spin only): `zeta = (rho_alpha - rho_beta) / (rho_alpha + rho_beta)`, in
  [-1, 1]: 0 for a closed shell or where there is no density, ±1 where only
  one spin is present.
- `shannon_entropy_density()` (`PropertyKind::ShannonEntropyDensity`, Total
  only): `-rho ln rho`, 0 where `rho <= 0`.
- `fisher_information_density()` (`PropertyKind::FisherInformationDensity`,
  Total only): `|grad rho|^2 / rho`, 0 where `rho <= 0`.

`Evaluator::evaluate` throws `std::invalid_argument` if a field is
requested for a component it is not defined for.

Validation (tests `density-fields` and `multiwfn-reference`):

- the values equal their formulas applied to the Evaluator's own value,
  gradient, Hessian, Laplacian and tau to 1e-13 (DORI 1e-12), for water,
  O2 (unrestricted) and acrolein;
- K agrees with the orbital form `-(1/2) sum_i psi_i lap psi_i`, evaluated
  with `Molecule::orbital_derivatives()`, to about 1e-13 of
  `|tau| + |lap rho| / 4`, and DORI with central finite differences of
  `|grad rho|^2 / rho^2` to about 1e-7 (the tests require 1e-7 and 1e-5);
- for a closed shell the total LOL equals the unpolarized formula and the
  spin polarization is 0; LOL and DORI lie in [0, 1] and `|zeta| <= 1`;
  2 and 4 threads give identical results;
- Multiwfn 2026.9.1 Cube files of water (LOL, IRI, DORI, K, V, H, Shannon,
  Fisher) and O2 (LOL, K, spin polarization) agree point by point within
  tolerances of 2e-5 to 1e-4, chosen per field (see
  [tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)); the largest
  difference is 41% of its tolerance.

## Orbital-energy fields: ALIE, LEA, LEAE, local electronegativity and hardness (1.13.0)

Five fields combine the SCF orbitals phi_i with their energies e_i
(hartree), with Multiwfn's definitions (its real-space function 18 and
user functions 27, -27, 28 and 29). With n_i = 2 for the orbitals of a
restricted closed-shell calculation and 1 for those of an unrestricted one:

- `average_local_ionization_energy()`: the ALIE of Sjoberg, Murray, Brinck
  and Politzer, `I = sum_occ n_i |e_i| |phi_i|^2 / rho`, over the occupied
  orbitals (rho is their density). Low values mark the electrons easiest to
  remove (electrophilic attack).
- `local_electron_affinity()`: `EA_L = -sum_vir e_i |phi_i|^2 / sum_vir
  |phi_i|^2` over all virtual orbitals (Brinck, Murray and Politzer, J. Mol.
  Model. 9, 342 (2003)). It depends on the high-lying virtuals, so it is
  meaningful only for basis sets without diffuse functions.
- `local_electron_attachment_energy()`: `E_att = sum_vir,e_i<0 n_i e_i
  |phi_i|^2 / rho`, over the virtual orbitals with negative energy only
  (Brinck, Carlqvist and Stenlid, J. Phys. Chem. A 120, 10023 (2016)). It
  is more robust than EA_L and admits diffuse functions, but it is 0
  everywhere when no virtual orbital has a negative energy (as for water at
  B3LYP/def2-TZVP). Negative values mark sites of nucleophilic attack; the
  original work uses the rho = 0.004 isosurface.
- `local_electronegativity()` and `local_hardness()`: `(I + EA_L) / 2` and
  `(I - EA_L) / 2`.

All five are defined for the total density only, use the orbitals whatever
the backend (with exact AO values, no screening), and are 0 where their
denominator vanishes, for example far from the molecule. They need an SCF
file whose occupied MOs reproduce its density matrices
(`ElectronicStructureInfo::ao_mo_eligible()`) and its orbital energies; the
last four also need virtual orbitals. Otherwise `evaluate()` throws
`std::invalid_argument`. The virtual orbitals are those of the FCHK above
the occupied ones (fewer than the basis functions if Gaussian removed
linear dependencies), so EA_L costs about as many orbital evaluations as
there are MOs. `average_local_ionization_energy(molecule, points)` (in
`<chem/kumo/surface.hpp>`) returns the same ALIE values, but throws where
the density vanishes.

Validation (tests `point-evaluation`, `multiwfn-reference`,
`molecular-surface`):

- the five fields equal their formulas evaluated from
  `Molecule::orbital_derivatives()` and the orbital energies to 1e-12 of
  the ALIE for water, O2 (unrestricted) and acrolein (LEAE also to 1e-12 of
  itself where |LEAE| >= 1e-3 hartree), and the ALIE equals
  `average_local_ionization_energy()`; LEAE is 0 for water, all fields are
  0 far from the molecule, and files without orbital energies or with a
  density the MOs do not reproduce are rejected;
- Multiwfn 2026.9.1 Cube files of acrolein (all five; one virtual orbital
  with negative energy) and O2 (ALIE, LEA, LEAE; two) agree point by point
  within 2e-5 relative. LEAE adds 1e-10 hartree absolute for its nodal
  planes and is held to 2e-5 relative alone where |LEAE| >= 1e-3 hartree
  (more than 10000 points per molecule; largest difference there 6.4e-6),
  and the electronegativity, which crosses zero, is compared on the scale
  (|I| + |EA_L|) / 2 (see
  [tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)). The largest
  difference is 54% of its tolerance.

## Electrostatic potential (1.9.0)

`density::total.electrostatic_potential()`
(`PropertyKind::ElectrostaticPotential`) is the molecular electrostatic
potential, in hartree/e:

    V(r) = sum_A Z_A / |r - R_A| - sum_(mu,nu) P_(mu nu) integral chi_mu(r') chi_nu(r') / |r - r'| dr'

with the nuclear charges Z_A of the FCHK (the `Nuclear charges` record when
present, so effective charges with ECPs) and the total AO density matrix P.
It is defined only for the total density; other components throw
`std::invalid_argument`. With the AO_MO backend the electronic term still
uses P, which then reproduces the orbital density (checked at load time).

At a grid point within 1e-8 bohr of a nucleus, that nucleus's own term,
singular there, is left out: the value is the potential that the nucleus
feels from the electrons and the other nuclei. Regular grids often contain
nuclear positions exactly: kumo's default grid for O2 does, for example.

The electronic term uses the McMurchie-Davidson scheme. The product of two
Cartesian Gaussian primitives is expanded in Hermite Gaussians at the
Gaussian product centre, and P is contracted with the expansion
coefficients once per evaluation. Each grid point then needs, for every
primitive pair, one Boys function F_m(T), the Hermite Coulomb integrals
R_tuv and a dot product. Pure shells are handled by transforming P to
Cartesian components. Shells up to l = 13 are supported.

Far field. The Evaluator processes the grid in clusters of 8 x 8 x 8
points and splits the error budget `potential_tolerance` (EvaluationOptions,
default 1e-8 hartree/e; 0 evaluates every pair exactly at every point)
between two approximations, each with a rigorous error bound:

- Source clusters (half of the budget). Primitive pairs are grouped by the
  2-bohr cell of their product centre and by exponent (a factor of 4 per
  group). The point-multipole moments of each group about its centre are
  computed once, and a group far enough from a cluster of points
  contributes through one multipole-to-local translation of total order 10.
  A group is used when the Taylor remainder of that translation and the
  neglected erfc terms of its Gaussians are each at most a quarter of the
  budget divided by the number of groups.
- Individual pairs (the other half). Each remaining pair has a bound of the
  error of its own Taylor expansion of order 6 about the cluster centre,
  tabulated at error levels 10^(-j/2) times the budget. In each cluster the
  pairs with the smallest bounds are expanded while the sum of their bounds
  fits in the budget; the rest are evaluated exactly at every point.

The bounds, in `src/electrostatic_potential.cpp`, split the potential of a
Hermite Gaussian into point-multipole (1/R) and erfc parts. The Taylor
remainders of the first use the bound m! / |x|^(m+1) on the derivatives of
order m of 1/|x|; the second is bounded through the integral representation
of erfc and Cramer's inequality for Hermite functions. The bounds are
conservative: on the test molecules the actual errors are at most 4e-10
with the default tolerance. Cluster boundaries depend only on the grid, so
the values do not depend on the block size or the number of threads.

Cost. Per grid point on one core: water (def2-TZVP, 556 primitive pairs)
about 23 microseconds exactly; HBC (def2-SVP, 177,315 pairs) about 3.7 ms
exactly, 0.50 ms with the default tolerance and 0.29 ms with 1e-6, averaged
over six sample clusters inside and around the molecule (Apple M, release
build; actual errors 4e-11 and 2e-8). The whole default grid of HBC
(1.5 million points) took 191 s on 8 cores with the default tolerance and
108 s with 1e-6; exactly it would take about 95 minutes on one core. The
error budget, not the arithmetic, limits the gain: with an absolute bound
of 1e-8 for a molecule with 234 electrons, the far field must be resolved to
about 1e-11 relative.

The Boys function uses a seventh-order Taylor expansion about a table with
spacing 0.05 for T < 40 and the closed form of F_0 with upward recursion
above, with downward recursion to the lower orders. It agrees with
quadrature of its definition to 2.3e-15 relative for m <= 32.

Validation (test `electrostatic-potential`, on water def2-TZVP, O2,
acrolein and HBC 6-31G(d,p), which has Cartesian d and SP shells):

- s-function densities on one and two centres agree with the erf closed
  form to 2e-15;
- 1000 bohr from the molecule, V agrees with the multipole expansion built
  from the charge, dipole and quadrupole of the analytic AO moments, with a
  residual below 1% of the quadrupole term (the octupole and higher terms);
- the Laplacian of the electronic term by fourth-order finite differences
  equals 4 pi rho within 8e-8 relative, the truncation error of the stencil
  (16 times smaller when the step is halved);
- with `potential_tolerance` 0 the Evaluator returns the same values as
  the exact pointwise engine, and with the default tolerance values within
  1e-8 of them (observed: at most 3e-10), for any block size, request and
  number of threads;
- clusters with tolerances 1e-4 and 1e-8, centred in and far from the
  molecule, stay within the tolerance (observed: at most 4e-6 and 4e-10).

The test `multiwfn-reference` also compares the ESP of water and O2 with
Multiwfn 2026.9.1 Cube files (its original analytic ESP code), point by
point: they agree within 2e-5 relative plus 1e-7 hartree/e, the rounding of
the reference files and of their grid coordinates (see
[tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)).

## Electrostatic potential on a molecular surface (1.10.0)

`<chem/kumo/surface.hpp>` computes the electrostatic potential on an
isosurface of the total electron density and the statistics of Politzer
and Murray that describe it (extrema, averages, variances, balance of
charges, internal charge separation, molecular polarity index):

```cpp
#include <chem/kumo.hpp>

const auto mol = chem::kumo::Molecule::from_fchk("acrolein.fchk.gz");
const auto surface = chem::kumo::electrostatic_potential_surface(mol);
const auto& s = surface.statistics;
std::cout << s.minimum << ' ' << s.maximum << ' ' << s.total_variance << ' '
          << s.balance << '\n';
for (std::size_t v : s.local_minima)  // vertices, by increasing value
    std::cout << surface.potential[v] << " at " << surface.mesh.vertices[v].transpose() << '\n';
```

`SurfaceOptions` sets the density isovalue (default 0.001 e/bohr^3, the
usual molecular surface), the grid spacing of the surface construction
(default 0.25 bohr) and the error bound of the potential at each vertex
(`potential_tolerance`, default 1e-8 hartree/e; 0 is exact).
`PotentialSurface` holds the mesh, the potential at its vertices (hartree/e)
and the statistics. The result does not depend on the number of OpenMP
threads.

Surface. `density_isosurface(molecule, isovalue, spacing)` returns the
isosurface alone as a `TriangleMesh`: vertex positions in bohr and
triangles oriented outwards, with `area()` and the enclosed `volume()`. The
total density is evaluated on an axis-aligned grid around the atoms, whose
margin grows from 4 bohr until the grid boundary lies outside the
isosurface. Each grid cell is split into six tetrahedra around its main
diagonal; the surface crosses every tetrahedron edge whose ends lie on
either side of the isovalue, and each tetrahedron contributes a triangle or
a quadrilateral (split along its shorter diagonal). Every vertex is then
moved onto the isosurface itself, by Newton steps along its edge on the
analytic density and its gradient, safeguarded by bisection, to 1e-10 bohr.
The mesh is closed and consistently oriented, so its volume follows from
the divergence theorem. Isovalues below 1e-8 e/bohr^3
(`minimum_density_isovalue`) are rejected: the AO cutoff changes the
density by about 1e-13 e/bohr^3. An isosurface that extends more than 16
bohr beyond the atoms is an error.

Potential. The vertices are evaluated in clusters (those on the edges of a
block of 8 x 8 x 8 grid cells) with the far-field approximations of the
Evaluator, within `potential_tolerance` at every vertex.

Statistics. `surface_statistics(mesh, values, polar_threshold)` works for
any function given at the vertices. It integrates the piecewise-linear
interpolation of the vertex values exactly over each triangle; the positive
region (V >= 0) and the negative one are bounded by the zero line of that
interpolation, which may cross triangles, and likewise for the averages of
|V - average| and |V| and for the nonpolar area (|V| <= `polar_threshold`,
default 10 kcal/mol). With A+ and A- the positive and negative areas:

- averages over the whole surface, A+ and A- (NaN for an empty region);
- variances sigma^2_+ and sigma^2_- about the positive and negative
  averages (0 for an empty region), and sigma^2_tot = sigma^2_+ + sigma^2_-;
- balance of charges nu = sigma^2_+ sigma^2_- / (sigma^2_tot)^2, in
  [0, 1/4] (NaN if sigma^2_tot = 0);
- internal charge separation Pi = (1/A) integral |V - average|, and the
  molecular polarity index MPI = (1/A) integral |V|, both in hartree/e;
- global extrema, the smallest and largest vertex values, and local
  extrema: vertices whose value is below (above) those of all neighbours
  and neighbours of neighbours, the criterion of Multiwfn. Where the
  potential is constant along a ring (linear molecules such as O2) the
  local extrema found on that ring depend on the mesh.

Accuracy and cost. The results converge at second order in the grid
spacing. On water, O2 and acrolein, area and volume at the default spacing
are within 0.15% of their values at 0.05 bohr, averages within 0.1% and
variances within 0.5% (O2) or 0.1%; global extrema within 2e-5 hartree/e.
Acrolein takes 1.3
s on one core; HBC def2-SVP (116,000 vertices) 28 s on 8 cores, of which
2.5 s for the surface.

Validation (tests `molecular-surface` and `multiwfn-surface-reference`):

- on a sphere, the vertices lie on it within 1e-9 bohr, the mesh is closed
  and oriented with Euler characteristic 2, area and volume converge to
  those of the sphere at second order, and the statistics of a linear
  function agree with their closed forms (averages R/2, variances R^2/12,
  nu = 1/4, Pi = MPI = R/2, nonpolar band 4 pi R c);
- on an octahedron the statistics of linear functions agree with
  hand-computed integrals to 1e-14, including regions cut through
  triangles;
- on water and O2 the vertices lie on rho = 0.001 within 1e-9 bohr, the
  potential agrees with the exact one within the tolerance (and to 1e-12
  with tolerance 0), and 2 and 4 threads give identical meshes and values;
- the statistics and the local extrema of water, O2 and acrolein agree with
  Multiwfn 2026.9.1 at a grid spacing of 0.05 bohr (see
  [tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)); the largest
  difference is 35% of its tolerance.

## Average local ionization energy on a molecular surface

The average local ionization energy (ALIE) of Sjoberg, Murray, Brinck and
Politzer is the orbital-energy average weighted by the orbital densities,

    I(r) = sum_i n_i |e_i| |phi_i(r)|^2 / sum_i n_i |phi_i(r)|^2,

over the occupied SCF orbitals phi_i with energies e_i (hartree) and
occupations n_i: 2 for the doubly occupied orbitals of a restricted
calculation, 1 otherwise (alpha and beta orbitals of an unrestricted one,
singly occupied orbitals of a restricted open-shell one). This is
Multiwfn's definition. Low values on the molecular surface mark the
electrons that are easiest to remove, so the surface minima indicate sites
of electrophilic attack.

```cpp
const auto surface = chem::kumo::local_ionization_energy_surface(mol);
const auto& s = surface.statistics;
for (std::size_t v : s.local_minima)  // by increasing value
    std::cout << surface.ionization_energy[v] * chem::kumo::hartree_in_electronvolt << " eV at "
              << surface.mesh.vertices[v].transpose() << " bohr\n";
```

`local_ionization_energy_surface(molecule, isovalue, spacing)` uses the
mesh of `density_isosurface()` (same defaults: 0.001 e/bohr^3 and 0.25
bohr) and returns an `IonizationEnergySurface`: the mesh, the energy at
each vertex (hartree) and its `SurfaceStatistics`. The ALIE is positive, so
the whole surface is the positive region: `average`, `total_variance` (the
variance about the average) and `charge_separation` (the average deviation)
are the meaningful statistics, with the extrema and the local extrema.
`average_local_ionization_energy(molecule, points)` evaluates it at any
points (columns in bohr).

The denominator is the density of the occupied orbitals. It equals the
FCHK total density to within the precision of the FCHK matrices, and only
when the occupied MOs reproduce the density matrices
(`ElectronicStructureInfo::ao_mo_eligible()`); otherwise, for example for a
post-SCF density, the functions throw `std::invalid_argument`, as they do
for a file without orbital energies, a nonfinite point or a point where the
density vanishes. Results do not depend on the number of OpenMP threads.

Accuracy and cost. On water, O2 and acrolein, the average at the default
spacing is within 1.5e-5 relative of Multiwfn at 0.05 bohr, the variance
within 0.6% (O2) and the global and local minima within 2.5e-4 hartree.
Local maxima lie on sharp ridges at the hydrogen tips; their number and
values depend on the mesh more than the minima do. Acrolein takes 0.3 s on
one core.

Validation (tests `molecular-surface` and `multiwfn-surface-reference`):

- at vertices and at points near and between the nuclei, the values agree
  with the orbital formula evaluated with `Molecule::orbital_derivatives()`
  to 1e-12 relative, for restricted (water) and unrestricted (O2) files;
- the values lie between the smallest and largest occupied |e_i|, and 2 and
  4 threads give identical results;
- files without orbital energies, or whose density matrix the occupied MOs
  do not reproduce, are rejected;
- the average, variance, global extrema and local minima of water, O2 and
  acrolein agree with Multiwfn 2026.9.1 at a grid spacing of 0.05 bohr (see
  [tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)); the largest
  difference is 35% of its tolerance.

## Any property on a molecular surface (1.13.0)

`map_on_density_isosurface(molecule, properties, isovalue, spacing,
options)` (in `<chem/kumo/surface.hpp>`) evaluates Evaluator properties at
the vertices of `density_isosurface()` with `Evaluator::evaluate(points)`
and returns a `MappedSurface`: the mesh, the requested properties (a
duplicate is mapped once), one value per vertex for each, and the
`SurfaceStatistics` of each (extrema, local extrema, averages, variances,
areas; the polar threshold is 0, so `nonpolar_area` carries no
information). Only one-component properties can be mapped (not `gradient()`
or `hessian()`); errors are those of `density_isosurface()` and
`Evaluator::evaluate()`, checked before the surface is built.

```cpp
using chem::kumo::density::total;
const auto surface = chem::kumo::map_on_density_isosurface(
    mol, {total.local_electron_attachment_energy()}, 0.004);
const auto& s = surface.statistics.front();  // LEAE minima: nucleophilic attack
```

The ALIE mapped this way equals `local_ionization_energy_surface()`
exactly; the electrostatic potential is evaluated as at any points and
agrees with `electrostatic_potential_surface()` within its tolerance (the
test `molecular-surface` finds at most 2e-8). Results do not depend on the
number of OpenMP threads.

`write_ply(path, mesh, fields, comments)` writes a mesh as an ASCII PLY
file: vertices (x, y, z in bohr) with one `double` property per field
(name, one value per vertex), then the triangles as vertex index lists.
Names must be letters, digits and `_`, distinct and different from x, y, z;
values must be finite; comments must not contain line breaks
(`std::invalid_argument` otherwise, before anything is written).

Validation (`multiwfn-surface-reference`): Multiwfn 2026.9.1's surface
analysis at 0.05 bohr, with kumo at its default 0.25 bohr, for the LEAE of
acrolein and O2 on rho = 0.004 (average within 3.0e-4 relative, variance
within 7.6e-3, global minimum within 4.6e-7 hartree, the six local minima of
acrolein within 1.8e-5 hartree and 0.07 bohr) and the LEA of acrolein on
rho = 0.001 (average 5.1e-5, variance 1.7e-3, local minima 1.9e-3 hartree,
which converge to Multiwfn's within 3e-5 at 0.05 bohr). The tolerances are
about three times these differences (see
[tests/REFERENCE_CUBES.md](../tests/REFERENCE_CUBES.md)). The two deepest
LEAE minima of acrolein lie on the beta carbon, above and below the
molecular plane, the site of nucleophilic (Michael) addition.

## Cube files (1.4.0)

`Field::write_cube(path, comment)` writes a scalar Gaussian Cube file in
bohr with six significant digits; non-finite values are rejected.
`read_cube(path)` (`<chem/kumo/cube.hpp>`) reads uncompressed or
gzip-compressed scalar Cube files (compression is detected from the gzip
signature) and returns a `CubeData`: a `Field` in bohr, whatever the units
of the header, and `CubeMetadata` with the two comment lines and the
original coordinate convention, which `write_cube(const CubeData&, path)`
keeps. Orbital Cube variants (negative atom counts, several datasets) are
not read.

## Auxiliary numerical APIs (1.4.0)

- `BasisEvaluator::overlap_matrix()` and `cartesian_moments()`: analytic AO
  overlap and Cartesian moments through second order about the origin, in
  FCHK AO order.
- `Molecule::density_matrix(spin)`: the AO density matrix of a component,
  as stored in (or derived from) the FCHK; the electron count is Tr(P S).
- `integrate_rectangular(field)`, `integrate_trapezoidal(field)` and
  `grid_cell_volume(field)` (`<chem/kumo/field_integrals.hpp>`): quadrature
  of a sampled field on its affine grid.
- `FieldSampler` (`<chem/kumo/field_sampler.hpp>`): trilinear interpolation
  of a Field at any point, returning `FieldSample{value, inside}` (zero
  outside the grid). The Field must outlive the sampler.
- `SubspaceProjector` (`<chem/kumo/subspace_projector.hpp>`): projector,
  orthogonal in the AO overlap metric, onto the span of a set of AO
  coefficient columns, with a relative rank tolerance.
- `WeightedOverlapAccumulator` and `weighted_overlap_matrix()`
  (`<chem/kumo/weighted_overlap.hpp>`): AO matrices weighted by a sampled
  field over a quadrature supplied by the caller.

## Errors and limitations

Malformed or unsupported FCHK records, invalid grids, invalid evaluation
options, properties requested for components they are not defined for, and
out-of-range MO indices, invalid surface options or meshes, and FCHK
files whose orbitals cannot give the average local ionization energy are
reported with exceptions. The library does not
infer absent MO coefficients or orbital energies. The reference data cover
O2, acrolein and water (Gaussian 09 Cubegen and Multiwfn Cube files), and
the HBC fixtures (6-31G(d,p) with Cartesian d and SP shells, def2-SVP,
def2-TZVP with f shells) cover larger bases; shells beyond f and other FCHK
variants (for example ECPs) have no reference cases yet. The mathematical
background is in [mathematics/chem-kumo-mathematics.tex](mathematics/chem-kumo-mathematics.tex).

The distributed FCHK fixtures use `.fchk.gz`; `Molecule::from_fchk` accepts both plain and gzip-compressed inputs.
