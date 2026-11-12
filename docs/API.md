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

## Errors and limitations

Malformed or unsupported FCHK records, invalid grids, and out-of-range MO
indices are reported with exceptions. The library does not infer absent MO
coefficients or orbital energies. Validation currently covers O2 and
acrolein reference files; higher-angular-momentum shells and other FCHK
variants require additional reference cases.

The distributed FCHK fixtures use `.fchk.gz`; `Molecule::from_fchk` accepts both plain and gzip-compressed inputs.
