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
- `densities_blocked(grid, components, block_size)` explicitly selects a
  block size (the default density path uses 128 points per block).
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

## Spatial derivative API (1.2.0)

`BasisEvaluator::evaluate_derivatives(point)` returns an `AODerivatives`
structure: `values` (N AO), `gradient` (N AO x 3) and `hessian`
(N AO x 6). The Hessian columns are **xx, yy, zz, xy, xz, yz**; the
coordinate axes and units are Cartesian bohr. AO ordering follows the FCHK.

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

`Molecule::density_derivatives_blocked(points, component, block_size=128)`
accepts a **3 x N** Eigen matrix with one Cartesian point per column,
and returns `std::vector<DensityDerivatives>` in the same order. A zero
block size and non-finite points are rejected. The method currently
materializes the result vector; it is not a streaming grid writer.

The blocked derivative path uses Eigen to apply the density matrix once
to a packed panel containing AO values and their three first derivatives.
Six AO second-derivative panels are used in linear reductions; no density
matrix product is needed for those panels. AO derivative evaluations
remain pointwise inside each block. There is no AO screening, sparse
matrix contraction, or new numerical dependency.

These are **spatial derivatives at fixed geometry and fixed wavefunction**,
not derivatives with respect to nuclear positions. The CLI has not been
extended to write derivative Cube files in this development branch.
The mathematical derivation is in
[chem-kumo-mathematics.tex](mathematics/chem-kumo-mathematics.tex).
Performance methodology and HBC measurements are in
[SPATIAL_DERIVATIVES.md](../benchmarks/SPATIAL_DERIVATIVES.md).

## Errors and limitations

Malformed or unsupported FCHK records, invalid grids, and out-of-range MO
indices are reported with exceptions. The library does not infer absent MO
coefficients or orbital energies. Validation currently covers O2 and
acrolein reference files; higher-angular-momentum shells and other FCHK
variants require additional reference cases.

The distributed FCHK fixtures use `.fchk.gz`; `Molecule::from_fchk` accepts both plain and gzip-compressed inputs.
