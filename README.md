# chem-kumo 1.1.0

**chem-kumo** is a lightweight C++20 library and command-line application for
reconstructing ground-state SCF electron densities and molecular orbitals from
Gaussian formatted checkpoint (`.fchk`) files. It supports restricted and
unrestricted calculations, with Eigen as its only numerical dependency.

## Requirements and build

- C++20 compiler; Eigen 3.4; Meson (>= 0.63) and Ninja for building.
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

## Command-line examples

```sh
kumo molecule.fchk                              # total SCF density
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
using namespace chem::kumo;

auto mol = Molecule::from_fchk("molecule.fchk");
auto grid = Grid::around(mol.atoms(), 5.0, 0.2); // bohr

mol.density(grid, Spin::Total).write_cube("density.cube");
auto spin_fields = mol.spin_densities(grid);
spin_fields.spin.write_cube("spin.cube");

// Orbital numbers are 1-based; Spin::Alpha / Spin::Beta select MO channels.
auto homo = mol.orbital(grid, 15, Spin::Alpha);
homo.write_orbital_cube("homo.cube", 15);
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

For a method-by-method guide and contracts, see [docs/API.md](docs/API.md).

## Validation and supported scope

The Meson suite currently registers **26 tests**: 13 established density,
input, and consumer tests; 10 orbital comparisons with Gaussian 09 Cubegen;
and 2 CLI listing checks. The 26 tests passed on macOS ARM64 (Apple Clang)
for the development package tested on 2026-10-13. The MO comparisons include
occupied and virtual alpha/beta orbitals of open-shell O2 and restricted
acrolein, including signed amplitudes and orbital norms. Existing density
references cover total, alpha, beta, and spin density. See
[tests/REFERENCE_CUBES.md](tests/REFERENCE_CUBES.md) for reference details.

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
