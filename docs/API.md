# chem-kumo C++ API guide

Include `<chem/kumo.hpp>` and link against `kumo_dep` in Meson. Public types
are in namespace `chem::kumo`. Examples use fully qualified names (such as
`chem::kumo::Molecule` and `chem::kumo::Spin`) rather than namespace-wide
`using` directives. All spatial coordinates and grid vectors use atomic
units (bohr).

## Reading an FCHK

`Molecule::from_fchk(path)` constructs the molecule, basis evaluator, SCF
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

## Errors and limitations

Malformed or unsupported FCHK records, invalid grids, and out-of-range MO
indices are reported with exceptions. The library does not infer absent MO
coefficients or orbital energies. Validation currently covers O2 and
acrolein reference files; higher-angular-momentum shells and other FCHK
variants require additional reference cases.
