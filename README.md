# chem-kumo 1.2.0

<p align="center">
  <img src="assets/chem-kumo-logo.png" alt="chem-kumo logo" width="280">
</p>

**chem-kumo** is a lightweight C++20 library and command-line application for
reconstructing ground-state SCF electron densities and molecular orbitals from
Gaussian formatted checkpoint (`.fchk` or gzip-compressed `.fchk.gz`) files. It supports restricted and
unrestricted calculations, with Eigen for numerical operations and zlib for gzip input.

## Requirements and build

- C++20 compiler; Eigen 3.4; zlib; Meson (>= 0.63) and Ninja for building.
- No Gaussian installation or external BLAS/LAPACK is needed to use KUMO.

```sh
./build-and-test.sh --clean
# Or manually:
meson setup build --buildtype=debugoptimized
meson compile -C build
meson test -C build --print-errorlogs
```

To install, use `meson install -C build` with a suitable Meson prefix.
The installable targets are the `kumo` executable, `libkumo`, and public headers.
Benchmarks are optional (`meson configure build -Dbenchmarks=true`).
For performance comparisons, select `-Dbuildtype=release` and rebuild.

## Command-line examples

```sh
kumo molecule.fchk                              # total SCF density
kumo molecule.fchk.gz                           # transparent gzip input
kumo o2.fchk --spin all                          # total, alpha, beta, spin difference
kumo molecule.fchk --mo homo                     # HOMO amplitude
kumo molecule.fchk --mo frontier                 # HOMO and LUMO
kumo molecule.fchk --mo 8,10:12,15               # explicit 1-based indices
kumo molecule.fchk --mo homo-2:lumo+2            # inclusive frontier-relative range
kumo o2.fchk --mo frontier --spin both           # alpha and beta frontiers
kumo molecule.fchk --list-mo                     # energies and occupations
kumo o2.fchk --list-mo --spin beta
kumo molecule.fchk --spacing 0.15 --margin 5 --output myfield
kumo --help
```

The default grid spacing is **0.20 bohr** and the margin around the atomic
bounding box is **5.0 bohr**. Both `--spacing` and `--margin` use bohr.
`--output` specifies a filename prefix, not a full Cube filename. Output
files are written to the current working directory unless the prefix includes
a directory. `--mo` accepts 1-based orbital numbers, comma-separated lists,
inclusive `A:B` ranges, `homo`, `lumo`, `frontier`, and offsets such as
`homo-2` and `lumo+2`. Indices refer to the chosen spin channel.

`--spin` selects a **density component** (`total`, `alpha`, `beta`,
`difference`, `all`) in density mode; with `--mo` or `--list-mo`, it selects
an **orbital channel** (`alpha`, `beta`, `both`). The default orbital channel
is alpha. Restricted calculations share one spatial MO set between alpha
and beta; unrestricted calculations have separate coefficients and frontier
indices. `--list-mo` displays occupation 2/0 for restricted spatial orbitals,
and 1/0 per spin channel for unrestricted orbitals.

Cube values for molecular orbitals are **signed amplitudes**
`psi_i(r) = sum_mu C[mu,i] chi_mu(r)`, **not** `|psi_i(r)|^2`. A global MO
phase change has no physical significance. Cube density values are in
electrons/bohr^3; orbital amplitudes are in bohr^(-3/2).

## C++ library API

```cpp
#include <chem/kumo.hpp>

int main()
{
    auto mol = chem::kumo::Molecule::from_fchk("molecule.fchk");
    auto grid = chem::kumo::Grid::around(mol.atoms(), 5.0, 0.2); // bohr

    mol.density(grid, chem::kumo::Spin::Total).write_cube("density.cube");
    auto spin_fields = mol.spin_densities(grid);
    spin_fields.spin.write_cube("spin.cube");

    // Orbital numbers are 1-based; Alpha/Beta select MO channels.
    auto homo = mol.orbital(grid, 15, chem::kumo::Spin::Alpha);
    homo.write_orbital_cube("homo.cube", 15);
    return 0;
}
```

The default density evaluator uses blocked Eigen contractions (128 points per
block). `density_pointwise()` is retained as a regression reference;
`densities_blocked()` allows an explicit block size. `orbitals()` evaluates
multiple MOs using shared AO batches. The API does not silently normalize
orbital fields or change their phase.

The grid origin, axis-step vectors, and atomic coordinates are in **bohr**.
`Grid::steps` stores the three axis vectors in its **columns**. Density
`Spin::Difference` means `rho_alpha - rho_beta`. The field index order is
`(i,j,k)` and a Cube file stores values on the specified regular grid.

### Spatial derivatives (since 1.2.0)

The C++ API evaluates analytic **spatial** derivatives at arbitrary points,
with nuclear geometry and FCHK wavefunction data held fixed. The
`BasisEvaluator` provides AO values, Cartesian gradients and all six unique
Hessian components. `Molecule::orbital_derivatives()` returns signed MO
amplitude, gradient and Hessian for a 1-based orbital index and alpha/beta
channel. `Molecule::density_derivatives()` returns density, gradient,
full Hessian, Laplacian and positive kinetic-energy density `tau` for any
density component. The blocked variant accepts a 3 x N point matrix:

```cpp
#include <chem/kumo.hpp>

int main()
{
    auto mol = chem::kumo::Molecule::from_fchk("molecule.fchk");
    const Eigen::Vector3d point(0.0, 0.0, 1.0); // bohr
    const auto d = mol.density_derivatives(point, chem::kumo::Spin::Total);
    const double laplacian = d.laplacian;
    const double kinetic_density = d.tau;
    (void)laplacian;
    (void)kinetic_density;
    return 0;
}
```

No derivatives with respect to nuclear displacements, SCF steps, or
response calculations are performed. The current `kumo` CLI still writes
density and MO Cube files; the new derivative methods are **C++ API only**.

For the mathematical definitions, see
[Mathematical Foundations](docs/mathematics/chem-kumo-mathematics.tex).
For numerical organization and performance measurements, see
[Spatial derivative benchmark](benchmarks/SPATIAL_DERIVATIVES.md).
For a method-by-method guide and contracts, see [docs/API.md](docs/API.md).

## Validation and supported scope

The Meson suite on the spatial-derivatives development branch registers
**30 tests**. The new tests cover AO, MO and density spatial derivatives
against centered finite differences; independent MO versus density-matrix
reconstruction; and blocked versus pointwise density derivatives. The 30
regression tests passed on macOS ARM64 (Apple Clang) on 2026-10-13 after
the fused Eigen optimization. Gaussian 09 Cubegen reference comparisons
for signed MO amplitudes and density components remain in the suite.
See [tests/REFERENCE_CUBES.md](tests/REFERENCE_CUBES.md) for reference details.

The reference tests account for the rounding of the **known** Gaussian grid
step (1/6 bohr) in Cube headers; this correction is specific to those
reference grids, not a general-purpose Cube-reading policy.

Validation is limited to the tested FCHK variants and basis shells; it is not
a guarantee of support for every Gaussian basis or every electronic-structure
method. The program evaluates the ground-state SCF data and MO coefficients
actually stored in the FCHK. Gaussian itself is not linked or called.

## Using KUMO as a Meson subproject

Put the project in `subprojects/chem-kumo` and add:

```meson
kumo_dep = dependency('chem-kumo', fallback: ['chem-kumo', 'kumo_dep'])
executable('myapp', 'main.cpp', dependencies: kumo_dep)
```

For a separate installation smoke test, run
`./tests/check-installed-consumer.sh`.

## Development

chem-kumo was developed by Matteo Tommasini in collaboration with OpenAI's
ChatGPT, through an iterative process of scientific design, C++20
implementation, numerical optimization, and validation against Gaussian
reference calculations. Scientific direction, testing, and final review were
carried out by Matteo Tommasini, with ChatGPT contributing to code development,
debugging, and documentation.

## License

MIT License; see [LICENSE](LICENSE). Eigen is an external dependency and is not
bundled. Reference Cube and FCHK files are provided for numerical regression.

### Compressed formatted checkpoints

KUMO reads gzip-compressed FCHK files directly with zlib; the CLI and
`chem::kumo::Molecule::from_fchk(path)` accept either plain `.fchk` or
compressed `.fchk.gz` without temporary decompressed files. Compression is
detected from the gzip header, not the filename. Corrupted gzip data and
invalid gzip checksums are reported as input errors. The FCHK parser and
all numerical algorithms are unchanged. zlib is a required build dependency.

To check the integrity of a committed compressed reference:

```sh
gzip -t tests/data/o2.fchk.gz
```

### Compressed reference data

The repository stores its five Gaussian formatted checkpoints as `.fchk.gz`
files: O2 and acrolein for numerical regression, and three HBC basis-set
variants in `tests/data/hbc/` for reproducible performance measurements.
KUMO reads them directly without temporary decompression. The `fchk-gzip`
regression test materializes plain-text copies in a temporary directory to
check exact numerical equivalence and rejects corrupted gzip trailers.
