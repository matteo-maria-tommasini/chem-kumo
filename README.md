# chem-kumo 1.12.0

<p align="center">
  <img src="assets/chem-kumo-logo.png" alt="chem-kumo logo" width="280">
</p>

**chem-kumo** is a lightweight C++20 library and command-line application for
reconstructing ground-state SCF electron densities and molecular orbitals from
Gaussian formatted checkpoint (`.fchk` or gzip-compressed `.fchk.gz`) files.
It supports restricted and unrestricted calculations, with Eigen for numerical
operations and zlib for gzip input. The C++ library also reads scalar Gaussian
Cube files, including gzip-compressed Cube input.

## Requirements and build

- C++20 compiler; Eigen 3.4 or 5 (both tested in CI); zlib; Meson (>= 0.63) and Ninja for building.
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

## Optional OpenMP parallelism (1.3.0)

OpenMP is optional and disabled by default. When enabled, independent grid
blocks are evaluated concurrently for electron densities, spatial density
derivatives, molecular orbitals, and blocked orbital derivatives. The
pointwise APIs remain serial. Since 1.7.0 this includes the Evaluator's
derivative and AO-MO path, about 4 times faster with 8 threads on HBC
def2-TZVP. Numerical kernels and point ordering are shared between serial
and parallel builds, and results do not depend on the thread count. Eigen internal threading is disabled in
OpenMP builds to avoid nested parallelism.

Use the build helper to select separate build directories and toolchains:

```sh
./build-and-test.sh --clean                 # serial, default compiler
./build-and-test.sh --clean --openmp --release  # OpenMP, release
```

On macOS, OpenMP builds use the supplied MacPorts Clang 22 native file
(`meson/native/macports-clang22.ini`) and require a compatible `libomp`.
On Linux, the default compiler (tested with GCC 14) is used. Manual Meson
configuration is also supported with `-Dopenmp=true` or `-Dopenmp=false`.
Requesting OpenMP without a suitable compiler/runtime fails at configuration.
Set `OMP_NUM_THREADS` to control thread count, and optionally use
`OMP_PROC_BIND`/`OMP_PLACES` for experiments with affinity; KUMO does not
impose a binding policy. See [OpenMP guide](docs/OPENMP.md).

## Command-line examples

```sh
kumo molecule.fchk                              # total SCF density
kumo molecule.fchk.gz                           # transparent gzip input
kumo o2.fchk --spin all                          # total, alpha, beta, spin difference
kumo molecule.fchk --property laplacian          # Laplacian of the total density
kumo molecule.fchk --backend ao-p                # force the density-matrix contraction
kumo molecule.fchk --property laplacian -v       # report progress on standard error
kumo molecule.fchk --property density,gradient,laplacian,tau
kumo o2.fchk --spin difference --property tau    # tau_alpha - tau_beta
kumo molecule.fchk --property elf                # electron localization function
kumo molecule.fchk --property rdg,sign-lambda2-rho  # NCI analysis
kumo molecule.fchk --property lol,iri,dori       # localization and interaction regions
kumo o2.fchk --spin difference --property spin-polarization
kumo molecule.fchk --property density,esp        # electrostatic potential
kumo molecule.fchk --property esp --esp-tolerance 1e-6  # faster, error <= 1e-6
kumo molecule.fchk --esp-surface                 # ESP statistics on rho = 0.001
kumo molecule.fchk --alie-surface                # ionization energy statistics
kumo molecule.fchk --property alie,lea,leae      # orbital-energy fields (Cube)
kumo molecule.fchk --point @1 --point 0,0,1.5 --property density,laplacian
kumo molecule.fchk --line @1:@2 --property density,lol   # profile -> PREFIX_line.dat
kumo molecule.fchk --plane @1:@2:@3 --property elf      # map -> PREFIX_plane.dat
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

`--property` (density mode only) selects what to write for each density
component, as a comma-separated list: `density` (default), `gradient` (the
norm |grad rho|, in e/bohr^4), `laplacian` (e/bohr^5) and `tau` (positive
kinetic energy density, hartree/bohr^3). Density files keep their names
(`PREFIX_density.cube`, `PREFIX_alpha.cube`, ...); derivative files are named
`PREFIX_<component>_<property>.cube`, for example
`PREFIX_density_laplacian.cube` (total density) or `PREFIX_spin_tau.cube`
(alpha minus beta).

Three derived properties are available too:
- `elf` is the electron localization function, in [0, 1], for `--spin total`,
  `alpha` or `beta`. The total ELF uses the spin-polarized formula.
- `rdg` is the reduced density gradient s, for the total density only.
  Values are capped at 100.
- `sign-lambda2-rho` is the total density times the sign of the middle
  eigenvalue of its Hessian, in e/bohr^3, for the total density only.

The `rdg` and `sign-lambda2-rho` files are the inputs of an NCI analysis:
an isosurface of s (for example 0.5), colored by sign(lambda2) rho. The
definitions are in [docs/API.md](docs/API.md).

More fields follow Multiwfn's definitions:
- `lol` is the localized orbital locator, in [0, 1], for `--spin total`,
  `alpha` or `beta`, like `elf`.
- `iri` (interaction region indicator, |grad rho| / rho^1.1, capped at 100)
  and `dori` (density overlap regions indicator, in [0, 1]) show covalent
  and weak interactions; total density only.
- `kinetic-k` is the kinetic energy density K = tau - lap rho / 4, for any
  `--spin`; `potential-energy` (V = lap rho / 4 - 2 tau) and
  `energy-density` (H = lap rho / 4 - tau) are for the total density.
- `spin-polarization` is (rho_alpha - rho_beta) / rho, in [-1, 1], for
  `--spin difference` only.
- `shannon` (-rho ln rho) and `fisher` (|grad rho|^2 / rho) are the Shannon
  entropy and Fisher information densities of the total density.
- `alie` (average local ionization energy), `lea` (local electron
  affinity), `leae` (local electron attachment energy),
  `local-electronegativity` and `local-hardness` weight the SCF orbitals by
  their energies (hartree); total density only. They need an SCF file whose
  occupied orbitals reproduce its density matrices, with orbital energies;
  `lea` uses all virtual orbitals and suits basis sets without diffuse
  functions, while `leae` uses only the virtual orbitals with negative
  energy and is 0 everywhere if there are none. See
  [docs/API.md](docs/API.md) for the formulas.

`esp` writes the molecular electrostatic potential (nuclear plus
electronic, in hartree/e) to `PREFIX_density_esp.cube`, for the total
density only. It is computed from Coulomb integrals over the AO density
matrix, with an error of at most 1e-8 hartree/e per point. Together with the density it gives the usual ESP map: an
isosurface of the density (for example 0.001 e/bohr^3) colored by the ESP.
At a grid point on a nucleus, which kumo's default grids often contain, the
nucleus's own singular term is left out, so the value is the potential at
that nucleus. `--esp-tolerance` sets the error bound at each point (default
1e-8 hartree/e; 0 is exact): far primitive pairs are approximated within
it. On the default grid, water (def2-TZVP) takes about 2 s on one core and
HBC (def2-SVP, 1.5 million points) about 3 minutes on 8 cores, or 2
minutes with `--esp-tolerance 1e-6`.

`--esp-surface` writes no Cube file: it computes the electrostatic potential
on the isosurface rho = 0.001 e/bohr^3 of the total density (`--isovalue`
changes it) and prints its statistics in hartree/e and kcal/mol: global and
local extrema (V_s,min and V_s,max, with positions in Angstrom), averages
over the whole surface and its positive and negative parts, variances
sigma^2_+, sigma^2_- and sigma^2_tot, the balance of charges nu, the
internal charge separation Pi, the molecular polarity index, and the
surface area, positive, negative and nonpolar (|V| <= 10 kcal/mol) areas
and enclosed volume. The surface is built on a grid with `--spacing`
(default 0.25 bohr here) and its vertices are placed on the analytic
isosurface; `--esp-tolerance` applies to the potential at the vertices.
Acrolein takes about 1 s; HBC (def2-SVP) about 30 s on 8 cores. The
results agree with Multiwfn's surface analysis (see
[docs/API.md](docs/API.md)).

`--alie-surface` does the same for the average local ionization energy, the
occupied orbital energies averaged with the orbital densities as weights,
in hartree and eV: global and local extrema (the minima mark the sites of
electrophilic attack), average, variance and average deviation, area and
volume. It needs the orbitals and orbital energies of an SCF calculation
whose occupied orbitals reproduce the FCHK density. Acrolein takes about
0.3 s.

`--point`, `--line` and `--plane` evaluate the `--property` and `--spin`
selection at chosen points instead of a Cube grid. Locations are `X,Y,Z` in
bohr or `@N`, the position of atom N (1-based). `--point` (repeatable)
prints a table on standard output. `--line A:B` writes `PREFIX_line.dat`,
evenly spaced points from A to B with both ends included (default spacing
0.01 bohr), extended by `--margin` beyond each end (default 0); its rows
give t, the distance from A towards B. `--plane A:B:C` writes
`PREFIX_plane.dat`, a rectangle in the plane of the three points covering
them with `--margin` on every side (defaults 0.10 bohr and 5.0 bohr); u
runs along A->B and v perpendicular to it towards C, and rows of constant u
are separated by blank lines, the layout gnuplot's `splot` expects. Each
table starts with `#` comment lines describing the geometry and every
column (quantity, units, component, contraction) and a `#` line with the
column names; the values have ten significant digits and the gradient is
written as its norm, as in Cube files. A plane map with gnuplot:

```gnuplot
set view map; splot 'molecule_plane.dat' using 1:2:6 with pm3d
```

`--backend` (density mode only) selects the density contraction: `auto`
(default) uses the occupied molecular orbitals (AO-MO) when they reproduce
the FCHK density matrices, which is checked when the file is loaded, and the
AO density matrices otherwise; `ao-p` always uses the density matrices, as
`kumo` did before 1.7.0; `ao-mo` requires the orbital contraction and fails
if the check did not pass. The check compares the total and spin density
matrices rebuilt from the occupied orbitals with those of the FCHK and
requires differences of at most 1e-8 relative to the total density matrix
for both; `--verbose` prints the two differences. The spin matrix is
measured against the total one because it is rebuilt as a difference of
two matrices of that size, so the rounding of the FCHK (about nine
significant digits) scales with the total matrix: measured against the
spin matrix itself, as in earlier versions, weakly spin-polarized states
such as broken-symmetry singlets could fail the check through rounding
alone. The two contractions agree to about 1e-9
relative, below the six significant digits written to Cube files (in tests,
at most the last printed digit of about 1% of the values changes), and
AO-MO is up to about 1.5 times faster for large basis sets. The second
comment line of each Cube file names the contraction used.

`-v` / `--verbose` reports on standard error what `kumo` is doing: the file
being read and a summary of the molecule (atoms, basis functions, SCF
reference and occupations), the result of the orbital check, the grid, the
density contraction and why it was chosen, each evaluation and writing step with its duration, and the
total time. When standard error is a terminal, evaluations also show a
progress bar with the elapsed and estimated remaining time. Standard output
and the Cube files are the same with and without `--verbose`.

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

### Typed Evaluator and electronic backends (1.4.0)

The C++ `chem::kumo::Evaluator` provides typed requests for total, alpha, beta,
and spin-difference electron density, gradients, Hessians, Laplacians and
positive kinetic-energy density (`tau`), the derived ELF, reduced density
gradient and sign(lambda2) rho (1.8.0), the electrostatic potential, as well
as orbital amplitudes. It
returns component arrays on the requested grid. This API is separate from the
existing `Molecule` convenience methods; the CLI uses it for all density
output.

`EvaluationOptions::backend` selects `Backend::AUTO` (default), `Backend::AO_P`
or `Backend::AO_MO`. `AUTO` prefers occupied-orbital AO-MO contractions when
orbital representation checks permit them; otherwise it uses the AO density
matrix (AO-P). An explicitly requested ineligible AO-MO backend is rejected.
`Results::backend_used()` reports the selected density backend. The AO-MO
performance advantage is workload-dependent, not a universal guarantee.

The Evaluator's `threads` option currently accepts `0` (runtime default) and
`1`; explicit values above one are rejected. For OpenMP-enabled builds, use
the runtime environment such as `OMP_NUM_THREADS` to control parallelism.
The grid evaluation uses bounded coordinate blocks controlled by `block_size`.

### AO screening, shell cutoff and performance (1.5.0, 1.6.0)

Grid densities and density derivatives skip, in each block of grid points,
the AOs that are negligible at every point of the block (value and
derivatives at most `default_screening_threshold` = 1e-12). The results
change only at rounding level; `EvaluationOptions::screening_threshold`
(and a parameter of `densities_blocked` and `density_derivatives_blocked`)
sets the threshold, and 0 gives the exact contraction. Together with batched
AO derivative evaluation, this makes 1.5.0 several times faster than 1.4.1
for large basis sets: on HBC def2-TZVP (1410 AOs), the total density is
about 10 times faster and density gradients and Hessians 4 to 6 times
faster. Since 1.6.0, each shell is also skipped at points beyond a
precomputed radius where a rigorous bound on its AOs and their first and
second derivatives is at most `BasisEvaluator::ao_cutoff` = 1e-14, which
makes densities a further 1.3 to 1.5 times and derivatives 1.1 to 1.3 times
faster on the same system. See [RELEASE_NOTES.md](RELEASE_NOTES.md) for the
measurements.

### Additional numerical APIs (1.4.0)

The library also provides analytic AO overlap integrals, spin-resolved
one-particle density matrices, metric subspace projectors, Cartesian AO moments,
trilinear sampling of scalar fields, and scalar-field integration. These
facilities have dedicated regression tests. See the installed public headers
in `include/chem/kumo/` for their exact types and contracts; not all are
available as command-line operations.

### Scalar Cube input (1.4.0)

`chem::kumo::read_cube(path)` reads uncompressed or gzip-compressed **scalar**
Gaussian Cube files. Compression is detected from the gzip signature rather
than the filename extension. The reader uses a common Cube parser, performs
incremental decompression without a temporary decompressed file, and rejects
invalid or truncated gzip streams. Its returned field uses bohr coordinates,
regardless of the coordinate convention in the Cube header.

The current Cube reader accepts nonnegative atom counts and one scalar dataset
(`NVAL` absent or equal to one); orbital Cube variants with negative atom counts
or multiple datasets are not supported by this reader. `write_cube()` writes
plain-text Cube output; gzip-compressed **output** is not part of this feature.
The `kumo` command-line application generates density and orbital Cube files
from FCHK input; it is not a general Cube inspection CLI.

## Validation and supported scope

The 1.12.0 release registers **81 Meson tests**, all passing in the
serial and OpenMP builds reported on macOS ARM64 on 2026-11-13, with both
Clang 22 and GCC 15, without compiler warnings. CI also runs the suite with
AddressSanitizer and UndefinedBehaviorSanitizer (except the three HBC cases
of the `slow` suite). Coverage
includes AO, MO and density spatial derivatives against finite differences;
independent MO versus AO-density-matrix reconstruction; backend selection;
grid validation at every evaluation entry point; AO screening against the
exact contraction; the AO shell cutoff bound; batched versus single-point AO
derivatives; command-line density and derivative Cube output, backend
selection and verbose mode; progress reports; blocked and pointwise
evaluation, including thread-count independence of the Evaluator and error
propagation from OpenMP workers; analytic
overlaps and moments; ELF, the reduced density gradient and
sign(lambda2) rho against their formulas, an independent orbital form of the
Pauli term and Multiwfn reference Cube files; the Boys function against
quadrature; the electrostatic potential against closed forms, the multipole
expansion of the analytic AO moments, the Poisson equation and Multiwfn,
and its far-field approximations against their error bounds; molecular
surfaces and their ESP and average local ionization energy statistics
against closed forms, the orbital formula and Multiwfn;
scalar-field operations; and
Cube and FCHK input. The Cube I/O tests exercise gzip equivalence and malformed
input. Gaussian 09 Cubegen reference comparisons for signed MO amplitudes and
density components remain in the suite.
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

## Compressed input and reference data

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

The test data include compressed Gaussian formatted checkpoints for O2,
acrolein, HBC basis-set variants and additional molecular validation cases.
The exact inventory is maintained under `tests/data/`. KUMO reads FCHK gzip
files directly without temporary decompression; the `fchk-gzip` regression
test compares compressed and plain-text inputs and rejects corrupted gzip
trailers. The scalar Cube I/O tests separately validate compressed Cube input.

## Development

chem-kumo was developed by Matteo Tommasini in collaboration with OpenAI's
ChatGPT, through an iterative process of scientific design, C++20
implementation, numerical optimization, and validation against Gaussian
reference calculations. Scientific direction, testing, and final review were
carried out by Matteo Tommasini, with ChatGPT contributing to code development,
debugging, and documentation.

The transition from version 1.4.0 to 1.4.1 resulted from a readability,
formatting and documentation audit conducted with Anthropic's Claude (via
Claude Code), which also implemented the resulting changes: consistent grid
validation, shared internal helpers, clearer error messages and test
diagnostics, and more robust CI. Every change was reviewed and approved by
Matteo Tommasini and verified against bit-identical numerical outputs.

The transition from version 1.4.1 to 1.5.0, a speed optimization of the
density and density-derivative evaluation, was likewise carried out with
Claude: profiling, a single density-matrix product per request, AO
screening, batched AO derivatives and AO-MO improvements. Matteo Tommasini
reviewed and approved every change; each optimization was checked against
the exact results (bit-identical, or rounding-level differences for AO
screening) and measured on the HBC reference systems.

Version 1.6.0, which adds density-derivative Cube output to the CLI and the
AO shell cutoff, was developed in the same way with Claude and reviewed and
approved by Matteo Tommasini, as was version 1.7.0 (CLI backend selection,
verbose mode with progress reports, and OpenMP parallelism of the Evaluator's
derivative path), version 1.8.0 (ELF, reduced density gradient and
sign(lambda2) rho, validated against Multiwfn), version 1.9.0 (the
molecular electrostatic potential with rigorous far-field error bounds),
version 1.9.1 (a code and documentation audit), version 1.10.0 (the
electrostatic potential on a molecular surface, validated against
Multiwfn), version 1.11.0 (the average local ionization energy on a
molecular surface, also validated against Multiwfn, and dynamic OpenMP
scheduling of the blocked grid loops) and version 1.12.0 (LOL, IRI, DORI,
the kinetic, potential and total energy densities, the spin polarization
and the Shannon and Fisher information densities, validated against
Multiwfn).

## License

MIT License; see [LICENSE](LICENSE). Eigen is an external dependency and is not
bundled. Reference Cube and FCHK files are provided for numerical regression.
