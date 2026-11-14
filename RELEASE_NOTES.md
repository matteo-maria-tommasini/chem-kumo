# chem-kumo 2.0.1 — release notes

A bug-fix release. Results and the public API are unchanged from 2.0.0;
apart from the version in the first line of Cube files, every CLI output
is byte-identical.

Fix:

- On machines with more OpenMP threads than blocks of points, the blocked
  evaluation loops allocated a workspace (the AO and product matrices of a
  block) in every thread, including threads that got no block. This wasted
  memory, which grows with the basis size, and made the `parallel-blocks`
  test fail with 176 threads ("made 176 for 63 blocks"). A thread now
  allocates its workspace when it gets its first block. The behaviour went
  back to 1.x; the test added in 1.13.1 exposed it.

Tests:

- `parallel-blocks` also checks 2 blocks and 0 blocks, which expose the
  old behaviour with any number of threads above the number of blocks. The
  whole OpenMP suite was run with `OMP_NUM_THREADS=176`.

# chem-kumo 2.0.0 — release notes

A new major version with a cleaner public API. libkumo had no external
users yet, so this was the time to remove the inconsistencies that had
accumulated over the 1.x series. The numerical results are those of 1.13.1:
apart from the version in the first line of Cube files, every CLI output is
byte-identical (checked on 529 output files, in serial and OpenMP builds),
and the `kumo` command line is unchanged.

Programs written for 1.x need changes; [docs/API.md](docs/API.md) opens
with a migration table from each 1.x function to its replacement.

Breaking changes:

- **Spin enumerations.** The four enumerations become two, in the new
  public header `<chem/kumo/spin.hpp>`: `DensityComponent {Total, Alpha,
  Beta, Difference}` selects a density, and `SpinChannel {Alpha, Beta}`
  selects the molecular orbitals (MOs) of one spin. `Spin` and
  `OrbitalSpin` are removed; `DensityComponent::Spin` is renamed
  `Difference` and `density::spin` is renamed `density::difference`.
  Asking for a "total" orbital no longer compiles.
- **One Hessian order.** Every six-component Hessian of the API lists the
  upper triangle row by row, `xx, xy, xz, yy, yz, zz`, as the Evaluator
  results always did. `AODerivatives::hessian`, the second-derivative
  panels of `BasisEvaluator::evaluate_derivatives_batch()` and
  `CartesianAOMoments::second` change (they were `xx, yy, zz, xy, xz,
  yz`); the values are the same, only their columns move.
- **Molecule holds the data; the Evaluator evaluates.** `Molecule` keeps
  loading and data: atoms, density matrices, MO coefficients
  (`mo_coefficients()` is now public), orbital energies and counts, and the
  orbital check. Its evaluation functions are removed from the public API:
  `density()`, `densities()`, `densities_blocked()`, `spin_densities()`,
  `density_pointwise()`, `density_derivatives()`,
  `density_derivatives_blocked()`, `orbital()`, `orbitals()`,
  `orbital_derivatives()` and `orbital_derivatives_blocked()`, with
  `SpinDensities`. `DensityDerivatives` and `OrbitalDerivatives` are no
  longer public. All evaluations go through the `Evaluator`, on grids or at
  points, with one set of options and the same validation.
- **Counts.** `Molecule::occupied_count()` returns `std::size_t`, like
  `orbital_count()`.

New:

- `orbital::alpha(i).gradient()` and `.hessian()` (and `beta`): the gradient
  and Hessian of an MO amplitude are Evaluator properties
  (`PropertyKind::OrbitalGradient`, `OrbitalHessian`), on grids or at
  points.
- `Evaluator::field(grid, property)` returns one one-component property as a
  `Field`, ready for `write_cube()`.
- On grids, the amplitudes of all requested orbitals of a spin channel are
  evaluated in one pass that shares the atomic-orbital values.

Tests and documentation:

- Tests and benchmarks reach the private kernels, and the simple pointwise
  reference implementations, through `tests/molecule_test_access.hpp`. A
  new test checks the Evaluator's orbital gradients and Hessians on a grid
  and at points against the pointwise reference (to 1e-12).
- `docs/API.md` is reorganized around the two main classes, with the
  migration table; README, `docs/OPENMP.md`, the benchmark notes and the
  mathematics document follow the new API.

# chem-kumo 1.13.1 — release notes

A maintenance release from a code, consistency and documentation audit.
Apart from the block-size fix below, results are unchanged: every CLI
output, error messages included, is byte-identical to 1.13.0 (checked on
529 output files, in serial and OpenMP builds). The public API is
unchanged.

Fix:

- `EvaluationOptions::block_size` is now checked once, when the
  `Evaluator` is built: it must be between 1 and the new constant
  `max_block_size` (65536, in `<chem/kumo/density.hpp>`), otherwise
  `std::invalid_argument` is thrown. Before, a block size above the limit
  was accepted by some evaluation paths and rejected by others, depending
  on the request.

Code organization (no change in results):

- The CLI is split into `src/kumo_main.cpp` (mode dispatch) and modules in
  `src/cli/`: option parsing, one table of the `--property` names with
  their components and methods, and one output module per mode. On a
  terminal, `--surface-map` now shows a progress bar like the other modes.
- `density.cpp` is split into `grid.cpp` (grids and `Field`),
  `molecule.cpp` (FCHK loading) and `density.cpp` (density kernels); the
  Gaussian-format Cube writer moves to `cube.cpp`.
- All block-parallel loops use two shared helpers; an allocation failure
  in a worker thread is now reported as an exception instead of
  terminating the program.
- The electrostatic potential on grids, on surface vertices and at
  arbitrary points shares one clustering implementation.
- `Evaluator::evaluate` and the isosurface extraction are split into named
  steps.
- The restricted open-shell code paths are kept on purpose, for future
  support of such files.

Documentation:

- The mathematics document is reorganized with one notation table and one
  meaning per symbol, and now covers the AO-MO contraction and its
  eligibility check, the automatic backend choice, the construction of
  grids, lines and planes, the clustering of the electrostatic potential
  and the values used where a field is undefined.
- `docs/OPENMP.md` lists every parallel loop and explains why results do
  not depend on the number of threads.
- `docs/API.md` and the README state that restricted open-shell (ROHF,
  ROKS) files cannot be loaded yet, since Gaussian writes no beta MO
  coefficients for them.
- Code comments throughout the library, the public headers and the CLI are
  rewritten for clarity, and about twenty that no longer matched the code
  are corrected (for example, the Boys-function Taylor expansion has degree
  six, seven terms).

Tests:

- New tests of the block-parallel helpers; the worker-error test now uses a
  throwing progress callback instead of an oversized block.

# chem-kumo 1.13.0 — release notes

The library and the `kumo` CLI evaluate density properties outside Cube
grids: at chosen points, along lines and on planes, and at the vertices of
the molecular surface. Four fields built from the virtual orbitals and
their energies join the average local ionization energy: the local electron
affinity, the local electron attachment energy, and the local
electronegativity and hardness. Any one-component property can be mapped on
the density isosurface, with its statistics and a PLY file of the colored
surface. The public API is extended in a source-compatible way; existing
outputs are unchanged.

Library changes:

- `Grid::line(start, end, margin, spacing)` and `Grid::plane(a, b, c,
  margin, spacing)` (in `<chem/kumo/density.hpp>`) build one- and
  two-dimensional grids, which every Evaluator property accepts.
- `Evaluator::evaluate(points, request)` evaluates any request at the
  columns of a 3 x N matrix (bohr), in column order. Values equal grid
  evaluation at the same points to rounding; the electrostatic potential
  groups the points by cubes of `potential_cluster_edge` (1.6 bohr) fixed
  in space, within `potential_tolerance`.
- New Total-only properties weighting the SCF orbitals by their energies,
  with Multiwfn's definitions: `average_local_ionization_energy()`,
  `local_electron_affinity()` (all virtual orbitals),
  `local_electron_attachment_energy()` (virtual orbitals with negative
  energy only; 0 everywhere if there are none),
  `local_electronegativity()` and `local_hardness()`. They are 0 where
  their denominator vanishes and need orbitals that reproduce the FCHK
  density matrices, with orbital energies (and virtual orbitals).
- `map_on_density_isosurface(molecule, properties, isovalue, spacing,
  options)` (in `<chem/kumo/surface.hpp>`) evaluates one-component
  properties at the vertices of the density isosurface and returns their
  values and `SurfaceStatistics`; `write_ply(path, mesh, fields,
  comments)` writes a mesh with vertex values as an ASCII PLY file.
- The surface ALIE (`average_local_ionization_energy(molecule, points)`)
  now shares the implementation of the new fields; its values are
  unchanged.

CLI changes:

- `--point X,Y,Z|@N` (repeatable) prints the `--property` and `--spin`
  selection at points, in bohr or at atom N; `--line A:B` writes
  `PREFIX_line.dat` (default spacing 0.01 bohr, margin 0) and `--plane
  A:B:C` writes `PREFIX_plane.dat` (default spacing 0.10 bohr, margin 5
  bohr, rows of constant u separated by blank lines for gnuplot). Tables
  start with `#` lines describing the geometry and every column.
- `--property` accepts `alie`, `lea`, `leae`, `local-electronegativity` and
  `local-hardness` (total density only).
- `--surface-map` prints the statistics of each `--property` on the density
  isosurface (extrema, averages, variances, areas, local extrema); `--ply`
  writes `PREFIX_surface.ply`, also with `--esp-surface` and
  `--alie-surface`, whose reports are unchanged. For the local electron
  attachment energy use `--isovalue 0.004`, the surface of the original
  work, which Multiwfn selects for it automatically.

Validation:

- New test `line-plane`: geometry and input checks, agreement with
  point-by-point evaluation, and LOL and ELF on a line and a plane through
  the points of the Multiwfn water Cube files (largest difference 0.24 of
  the tolerance).
- New test `point-evaluation`: every property at points against grid
  evaluation (water, O2, acrolein; both backends; at most 5e-15 relative,
  2e-10 for the electrostatic potential), the new fields against their
  orbital formulas (1e-12), identical results with 1, 2 and 4 threads.
- `multiwfn-reference`: eight new Multiwfn 2026.9.1 Cube files (acrolein:
  all five fields; O2: ALIE, LEA, LEAE) agree within 2e-5 relative; LEAE
  must meet that tolerance alone on more than 10000 points per molecule
  where |LEAE| >= 1e-3 hartree.
- `multiwfn-surface-reference`: LEAE on rho = 0.004 (acrolein, O2) and LEA
  of acrolein against Multiwfn's surface analysis, with kumo at 0.25 bohr
  and Multiwfn at 0.05 bohr; the six LEAE minima of acrolein (the deepest
  on the beta carbon) agree within 1.8e-5 hartree. Largest difference 0.36
  of its tolerance.
- `molecular-surface`: the mapped ALIE equals the ALIE surface analysis
  exactly; PLY files are read back; invalid input is rejected.
- 83/83 Meson tests passed on macOS ARM64 on 2026-11-13 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings. CI builds and tests on ubuntu-latest
  and macos-latest, and runs the suite under AddressSanitizer and
  UndefinedBehaviorSanitizer without the `slow` suite.

Documentation: the README, `docs/API.md`, the mathematics document and
`tests/REFERENCE_CUBES.md` describe the new features, their definitions,
edge values and validation.

---

# chem-kumo 1.12.0 — release notes

The library and the `kumo` CLI compute nine more real-space fields derived
from the electron density, with Multiwfn's definitions: the localized
orbital locator (LOL), the interaction region indicator (IRI), the density
overlap regions indicator (DORI), the kinetic, potential and total energy
densities K, V and H, the spin polarization and the Shannon and Fisher
information densities. The public API is extended in a source-compatible
way; existing outputs are unchanged.

Library changes (`<chem/kumo/evaluator.hpp>`):

- New `PropertyKind` values with `DensityDescriptor` helpers: `Lol`
  (`lol()`), `InteractionRegionIndicator`, `DensityOverlapRegionsIndicator`,
  `HamiltonianKineticEnergy`, `PotentialEnergyDensity`, `EnergyDensity`,
  `SpinPolarization`, `ShannonEntropyDensity`, `FisherInformationDensity`.
- LOL = t / (1 + t), t = D0 / tau, for Total, Alpha or Beta; the total LOL
  sums tau and D0 over both channels, as the total ELF does.
- IRI = |grad rho| / rho^1.1, capped at `interaction_region_indicator_cap`
  (100); DORI = theta / (1 + theta), computed from q = grad rho / rho without
  dependence on the scale of rho.
- K = tau - lap rho / 4 (any component), V = lap rho / 4 - 2 tau and
  H = lap rho / 4 - tau = -K (Total).
- zeta = rho_spin / rho (Spin component), -rho ln rho and |grad rho|^2 / rho
  (Total).
- Edge values are documented (no density, critical points, tau = 0), and
  `Evaluator::evaluate` throws `std::invalid_argument` for a component a
  field is not defined for. Results do not depend on the number of OpenMP
  threads.

CLI changes:

- `--property` accepts `lol`, `iri`, `dori`, `kinetic-k`,
  `potential-energy`, `energy-density`, `spin-polarization`, `shannon` and
  `fisher`, each written to its own Cube file with a descriptive comment.
- Combinations a field is not defined for (for example `iri` with
  `--spin alpha`, `lol` with `--spin difference`, `spin-polarization`
  without it) are rejected, and an unknown property name lists the valid
  ones.

Validation:

- The new test `density-fields` checks each field against its formula
  applied to the Evaluator's own value, gradient, Hessian, Laplacian and
  tau to 1e-13 (DORI 1e-12) for water, O2 and acrolein; K against the
  orbital form -(1/2) sum_i psi_i lap psi_i to about 1e-13; DORI against
  central finite differences to about 1e-7; closed-shell limits, ranges,
  thread independence, far-field values and rejected components.
- The test `multiwfn-reference` compares new Multiwfn 2026.9.1 Cube files
  of water (LOL, IRI, DORI, K, V, H, Shannon, Fisher) and O2 (LOL, K, spin
  polarization) point by point, with tolerances of 2e-5 to 1e-4 chosen per
  field and documented in `tests/REFERENCE_CUBES.md`; the largest
  difference is 41% of its tolerance. No existing tolerance was changed.
- `cli-end-to-end` checks the ranges of LOL, DORI, IRI and zeta, H = -K,
  and the rejected combinations.
- 81/81 Meson tests passed on macOS ARM64 on 2026-11-13 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings. CI builds and tests on ubuntu-latest
  and macos-latest, and runs the suite under AddressSanitizer and
  UndefinedBehaviorSanitizer without the `slow` suite.

Documentation: `docs/API.md`, the README and the mathematics document
describe the new fields and their edge values.

---

# chem-kumo 1.11.0 — release notes

The library and the `kumo` CLI compute the average local ionization energy
(ALIE) of Sjoberg, Murray, Brinck and Politzer, at any points and on the
molecular surface, as a companion to the electrostatic potential surface
analysis of 1.10.0. The blocked grid loops of OpenMP builds are scheduled
dynamically and run faster. The public API is extended in a
source-compatible way; existing outputs are unchanged.

Library changes:

- `average_local_ionization_energy(molecule, points)` (in
  `<chem/kumo/surface.hpp>`): the ALIE in hartree at the columns of
  `points`, sum_i n_i |e_i| |phi_i|^2 / sum_i n_i |phi_i|^2 over the
  occupied SCF orbitals, with occupations 2 for the doubly occupied
  orbitals of a restricted calculation and 1 otherwise; this is Multiwfn's
  definition.
- `local_ionization_energy_surface(molecule, isovalue, spacing)` returns an
  `IonizationEnergySurface`: the mesh of `density_isosurface()` (defaults
  0.001 e/bohr^3 and 0.25 bohr), the ALIE at its vertices and its
  `SurfaceStatistics` (extrema and local extrema, average, variance,
  average deviation). `hartree_in_electronvolt` converts the energies.
- Both throw `std::invalid_argument` for a file without orbital energies or
  whose occupied MOs do not reproduce the density matrices (for example a
  post-SCF density), a nonfinite point, or a point where the density
  vanishes. Results do not depend on the number of OpenMP threads.

CLI changes:

- `kumo FILE --alie-surface` prints the ALIE analysis in hartree and eV,
  with the local minima (sites of electrophilic attack) and maxima and
  their positions in Angstrom. `--isovalue` and `--spacing` apply as for
  `--esp-surface`; options that only concern Cube output, and
  `--esp-tolerance`, are rejected. The `--esp-surface` report is unchanged.

Performance:

- OpenMP builds assign the blocks of the blocked density, density
  derivative, AO-MO and orbital loops with `schedule(dynamic)` instead of
  `schedule(static)`: with AO screening, blocks far from the molecule cost
  much less than those near it, and a static assignment left threads idle.
  On the Apple M3 (4P+4E), HBC with 103,823 points, the blocked density
  speedup at 4 / 8 threads rose from 2.8x / 4.0x to 3.7x / 5.0x
  (def2-SVP) and from 2.5x / 3.7x to 3.5x / 4.8x (def2-TZVP). Results are
  bit-identical across schedules and thread counts.

Validation:

- The test `molecular-surface` checks the ALIE against the orbital formula
  evaluated with `Molecule::orbital_derivatives()` to 1e-12 relative, for
  restricted (water) and unrestricted (O2) files, its bounds, thread
  independence and the rejected inputs.
- The test `multiwfn-surface-reference` compares water, O2 and acrolein
  with new Multiwfn 2026.9.1 ALIE references at 0.05 bohr (kumo at its
  default 0.25 bohr): average within 1.5e-5 relative, variance within 6e-3
  (O2), global extrema and local minima within 2.5e-4 hartree. The
  tolerances are about three times these discretization differences,
  documented in `tests/REFERENCE_CUBES.md`; the largest difference is 35%
  of its tolerance.
- 80/80 Meson tests passed on macOS ARM64 on 2026-11-13 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings. CI builds and tests on ubuntu-latest
  and macos-latest, and runs the suite under AddressSanitizer and
  UndefinedBehaviorSanitizer without the `slow` suite.

Documentation: `docs/API.md`, the README and the mathematics document
describe the ALIE; `docs/OPENMP.md` and `benchmarks/PARALLEL_DENSITY.md`
have new release-build speedups of the blocked density and of the
`Evaluator` point path on the Apple M3.

---

# chem-kumo 1.10.0 — release notes

The library and the `kumo` CLI analyse the electrostatic potential (ESP) on
a molecular surface, an isosurface of the total electron density, with the
statistics of Politzer and Murray. The public API gains a header and is
otherwise source-compatible; `Molecule` has a new friend declaration, so
rebuild code that links against libkumo. Existing outputs are unchanged.

Library changes:

- New header `<chem/kumo/surface.hpp>` (included by `<chem/kumo.hpp>`):
  `electrostatic_potential_surface(molecule, SurfaceOptions)` returns a
  `PotentialSurface` with the mesh, the potential at its vertices
  (hartree/e) and its `SurfaceStatistics`. `SurfaceOptions` sets the
  density isovalue (default 0.001 e/bohr^3), the grid spacing (default
  0.25 bohr) and the potential tolerance (default 1e-8 hartree/e).
- `density_isosurface(molecule, isovalue, spacing)` returns the isosurface
  alone as a `TriangleMesh` (vertices in bohr, outward-oriented triangles,
  `area()` and `volume()`). It uses marching tetrahedra (six per grid cell)
  and moves every vertex onto the analytic isosurface by safeguarded Newton
  steps, to 1e-10 bohr. The mesh is closed and consistently oriented.
  Isovalues below 1e-8 e/bohr^3 (`minimum_density_isovalue`) are rejected,
  as is a surface extending more than 16 bohr beyond the atoms.
- `surface_statistics(mesh, values, polar_threshold)` works for any
  function given at the vertices and integrates its piecewise-linear
  interpolation exactly, cutting triangles along V = 0 and |V| =
  threshold: global and local extrema (Multiwfn's criterion), averages over
  the surface and its positive and negative parts, the variances
  sigma^2_+, sigma^2_- and sigma^2_tot, the balance of charges nu, the
  internal charge separation Pi, the molecular polarity index MPI, and the
  total, positive, negative and nonpolar areas.
- The potential at the vertices uses the far-field approximations of the
  Evaluator, within the tolerance at every vertex; results do not depend
  on the number of OpenMP threads.

CLI changes:

- `kumo FILE --esp-surface` prints the analysis in hartree/e and kcal/mol,
  with areas in bohr^2 and Angstrom^2, the enclosed volume and the local
  extrema with their positions in Angstrom. It writes no Cube file.
  `--isovalue` sets the density isovalue, and `--spacing` (default 0.25
  bohr with `--esp-surface`) and `--esp-tolerance` apply as well.
  Options that only concern Cube output are rejected with `--esp-surface`.

Validation:

- New test `molecular-surface`: on a sphere the mesh is closed with Euler
  characteristic 2, its vertices lie on the sphere within 1e-9 bohr, area
  and volume converge at second order and the statistics of a linear
  function agree with their closed forms; on an octahedron the statistics
  agree with hand-computed integrals to 1e-14; on water and O2 the
  vertices lie on rho = 0.001 within 1e-9 bohr, the potential agrees with
  the exact one and 2 and 4 threads give identical results.
- New test `multiwfn-surface-reference`: water, O2 and acrolein against
  Multiwfn 2026.9.1 at 0.05 bohr (kumo at its default 0.25 bohr). The
  tolerances are about three times the observed discretization
  differences, documented in `tests/REFERENCE_CUBES.md`; the largest
  difference is 35% of its tolerance.
- 80/80 Meson tests passed on macOS ARM64 on 2026-11-13 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings. CI builds and tests on ubuntu-latest
  and macos-latest, and runs the suite under AddressSanitizer and
  UndefinedBehaviorSanitizer without the `slow` suite.

Documentation: `docs/API.md`, `docs/OPENMP.md`, the README and the
mathematics document describe the surface construction, the statistics
and their accuracy.

---

# chem-kumo 1.9.1 — release notes

A maintenance release after a code and documentation audit. There are no
API or CLI changes, and results are identical to 1.9.0: every Cube file of
a 90-file regression set (densities, derivatives, orbitals, derived
properties and the electrostatic potential, serial and OpenMP) is
byte-identical. Public headers changed only in comments and a private
member, so rebuild code that links against libkumo.

Fixes:

- `Field::write_cube`, which the CLI uses, now throws
  `std::invalid_argument` for NaN or infinite samples instead of writing a
  Cube file that `read_cube()` and other readers reject, as
  `write_cube(CubeData)` already did.

Documentation:

- `docs/API.md` has new sections on the Evaluator (requests, results,
  options, lifetime, errors, determinism), Cube input and output and the
  auxiliary numerical APIs, and an updated list of errors and limitations.
- `docs/OPENMP.md` describes the derived properties and the electrostatic
  potential pass.
- The mathematics document (`docs/mathematics/`) covers AO screening and
  the shell cutoff bound, ELF, the reduced gradient, sign(lambda2) rho and
  the electrostatic potential (Hermite expansion, Boys function, far-field
  bounds and error budget).
- Public headers document the types, members and functions that had no
  comment.

Internal changes:

- Shared helpers replace duplicated code: one parallel loop over blocks,
  the normalized AO blocks of a basis, the density-derivative contraction
  kernel, grid point indexing, and smaller pieces of the orbital, backend
  choice and Cube I/O code. The electrostatic potential far-field code is
  split into planning and evaluation helpers.

Validation:

- Every test failure now reports what was expected, through
  `tests/test_support.hpp`.
- The parallel density, derivative and orbital tests require results
  identical to the single-thread run, bit for bit, at 2 and 4 threads with
  the same block size.
- CI compiles the benchmarks (`-Dbenchmarks=true`), so an API change cannot
  leave them broken, and the HBC electrostatic potential test has a 300 s
  timeout for slower hosted runners.
- 78/78 Meson tests passed on macOS ARM64 on 2026-11-13 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings. CI builds and tests on ubuntu-latest
  and macos-latest, and runs the suite under AddressSanitizer and
  UndefinedBehaviorSanitizer without the `slow` suite.

---

# chem-kumo 1.9.0 — release notes

The Evaluator and the `kumo` CLI compute the molecular electrostatic
potential (ESP), with exact Coulomb integrals near each point and far-field
approximations within an explicit error bound. The public API is extended
in a source-compatible way; rebuild code that links against libkumo
(`PropertyKind` has a new value and `EvaluationOptions` a new member).
Existing outputs are unchanged.

Library changes:

- `density::total.electrostatic_potential()`
  (`PropertyKind::ElectrostaticPotential`): the nuclear plus electronic
  potential in hartree/e, from the FCHK nuclear charges (the `Nuclear
  charges` record when present, so effective charges with ECPs) and the
  total AO density matrix. Total density only; other components throw
  `std::invalid_argument`.
- The electronic term uses McMurchie-Davidson Hermite expansions and a new
  Boys function (Taylor table below T = 40, closed form above; agrees with
  quadrature to 2.3e-15). Shells up to l = 13 are supported.
- At a grid point within 1e-8 bohr of a nucleus, that nucleus's own
  singular term is left out, so the value is the potential at the nucleus.
  Regular grids, including kumo's default grids, often contain nuclei.
- `EvaluationOptions::potential_tolerance` (default 1e-8 hartree/e) bounds
  the error at each point. The grid is processed in clusters of 8 x 8 x 8
  points. Groups of primitive pairs far from a cluster enter through
  multipole-to-local translations of order 10, other far pairs through
  Taylor expansions of order 6 about the cluster centre, and the rest
  exactly. The sum of rigorous error bounds stays within the tolerance; 0
  evaluates every pair exactly. Values do not depend on the block size or
  the number of threads.

Command-line changes:

- `--property esp` writes `PREFIX_density_esp.cube` (total density only).
  The Cube comment states the error bound, or "exact Coulomb integrals".
- `--esp-tolerance HARTREE` sets the bound (default 1e-8; 0 is exact). It
  is accepted only with `--property esp`.

Performance:

- HBC def2-SVP (177,315 primitive pairs) on its default grid (1.5 million
  points) takes 191 s on 8 cores with the default tolerance, or 108 s with
  1e-6. Exactly it would take about 95 minutes on one core. Water
  def2-TZVP takes about 2 s on one core.

Validation:

- New test `boys-function` against quadrature, closed forms and the
  recurrence.
- New test `electrostatic-potential`:
  - s-function densities agree with the erf closed form to 2e-15;
  - 1000 bohr away, V matches the multipole expansion from the analytic AO
    moments within 1% of the quadrupole term;
  - the Poisson equation holds within the 8e-8 truncation error of the
    finite-difference stencil;
  - far-field approximations stay within their tolerance, with observed
    errors at most 4e-10 at 1e-8;
  - HBC 6-31G(d,p) covers Cartesian d and SP shells (`slow` suite).
- `multiwfn-reference` compares the ESP of water and O2 with Multiwfn
  2026.9.1 (its original analytic code, `iESPcode=1`), within 2e-5 relative
  plus 1e-7 hartree/e, the rounding of the reference Cube files.
- 78/78 Meson tests passed on macOS ARM64 on 2026-11-12 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings. CI builds and tests on ubuntu-latest
  and macos-latest, and runs the suite under AddressSanitizer and
  UndefinedBehaviorSanitizer without the `slow` suite.

---

# chem-kumo 1.8.0 — release notes

The Evaluator and the `kumo` CLI compute three derived density properties:
the electron localization function (ELF), the reduced density gradient (s)
and sign(lambda2) rho, the inputs of an NCI analysis. The public API is
extended in a source-compatible way; rebuild code that links against
libkumo (`PropertyKind` has new values). Existing outputs are unchanged.

Library changes:

- `density::total.elf()`, `density::alpha.elf()`, `density::beta.elf()`
  (`PropertyKind::Elf`): the Becke–Edgecombe ELF, in [0, 1]. The total ELF
  uses the spin-polarized form, summing the Pauli kinetic energy density
  and its uniform-gas value over the alpha and beta channels; for a closed
  shell this equals the classic formula. ELF is 0 where there is no
  density, and is not defined for the spin difference.
- `density::total.reduced_gradient()` (`PropertyKind::ReducedGradient`):
  `s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3))`, capped at the new
  constant `reduced_gradient_cap` (100), also where rho <= 0, so Cube files
  stay finite.
- `density::total.sign_lambda2_density()`
  (`PropertyKind::SignLambda2Density`): rho times the sign of the middle
  eigenvalue of the density Hessian.
- The reduced gradient and sign(lambda2) rho are defined only for the total
  density; `Evaluator::evaluate` throws `std::invalid_argument` for an
  undefined component.
- The properties use the existing backends, AO screening and OpenMP
  scheduling; results do not depend on the number of threads. ELF and s
  cost about as much as tau, and sign(lambda2) rho about as much as the
  Laplacian.

Command-line changes:

- `--property elf` (with `--spin total`, `alpha` or `beta`), `--property rdg`
  and `--property sign-lambda2-rho` (total density only) write
  `PREFIX_density_elf.cube`, `PREFIX_density_rdg.cube`,
  `PREFIX_density_sign-lambda2-rho.cube` and, for spin channels,
  `PREFIX_alpha_elf.cube` and so on. Unsupported spin choices fail with a
  clear message.

Validation:

- The derived values equal their formulas applied to the Evaluator's
  density derivatives to about 1e-12; the Pauli term agrees with the
  cancellation-free orbital (Lagrange) form to better than 1e-9 relative;
  default AO screening changes ELF by at most 2e-8 and s by at most 5e-10
  relative on the O2 and acrolein grids.
- New test `multiwfn-reference`: ELF, s and sign(lambda2) rho of water and
  the ELF of O2 (open shell) agree with Multiwfn 2026.9.1 Cube files within
  1.2e-5 relative, at the six significant digits of the references
  (`tests/data/multiwfn/`, provenance in `tests/REFERENCE_CUBES.md`).
- On water, the O-H bond critical point has s near 0, sign(lambda2) rho =
  -0.361 e/bohr^3 and ELF 0.98.
- CI now also runs the test suite with AddressSanitizer and
  UndefinedBehaviorSanitizer on ubuntu-latest, and fails on compiler
  warnings (`-Dwerror=true`). The HBC cases of `density-screening` and
  `ao-shell-cutoff` are separate tests in a `slow` suite, which the
  sanitizer job skips.
- 75/75 Meson tests passed on macOS ARM64 on 2026-11-12 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings; CI builds and tests on ubuntu-latest
  and macos-latest.

---

# chem-kumo 1.7.1 — release notes

A maintenance release: release builds with GCC 14 and 15 are now
warning-free. There are no API, CLI or numerical changes; results are
identical to 1.7.0.

- Eigen is now a system dependency (`include_type: 'system'`, i.e.
  `-isystem`). GCC 14 (Debian) and GCC 15 reported many
  `-Wmaybe-uninitialized` and `-Wclass-memaccess` warnings inside Eigen 3.4
  headers, known false positives of the optimizer; warnings raised inside
  Eigen's headers are no longer reported, while warnings in this project's
  code still are. The generated pkg-config file is unchanged.
- `tests/cube_io_test.cpp` kept a lazy Eigen expression (`const auto`) that
  referenced a temporary `Eigen::Vector3d` after the statement; the point
  is now stored as a concrete `Eigen::Vector3d`. This affected only the
  test.
- The loop counter in the internal `factorial()` of the basis evaluator is
  now `long long`, so it cannot overflow (GCC 15
  `-Waggressive-loop-optimizations`); the multiplication order and the
  results are unchanged.

Validation:

- 71/71 Meson tests passed on macOS ARM64 on 2026-11-12 with MacPorts
  Clang 22 (serial and OpenMP) and GCC 15 (serial and OpenMP release
  builds), all four without warnings; CI builds and tests on ubuntu-latest
  and macos-latest.

---

# chem-kumo 1.7.0 — release notes

The command-line tool chooses the faster density contraction by default and
can report its progress, and the Evaluator's derivative and AO-MO path now
runs on all OpenMP threads. The public API is extended in a
source-compatible way; rebuild code that links against libkumo (the layout
of `EvaluationOptions` changed).

Command-line changes:

- `--backend auto|ao-p|ao-mo` selects the density contraction, and the
  default is now `auto`, as in the Evaluator: the occupied-orbital
  contraction (AO-MO) is used when the FCHK orbitals reproduce its density
  matrices (checked when the file is loaded), and the AO density matrices
  (AO-P) otherwise. **Behaviour change:** before 1.7.0 the CLI always used
  AO-P. The two contractions agree to about 1e-9 relative, below the six
  significant digits written to Cube files; on O2 and acrolein at most the
  last printed digit of about 1% of the values changes. `--backend ao-p`
  reproduces the Cube files of 1.6.0 digit for digit; `--backend ao-mo`
  fails clearly for an ineligible FCHK. The second comment line of every
  density Cube file names the contraction used.
- `-v` / `--verbose` reports on standard error the file being read, a
  summary of the molecule (atoms, basis functions, SCF reference and
  occupations), the grid, the density contraction and why it was chosen,
  each evaluation and writing step with its duration, and the total time.
  On a terminal, evaluations show a progress bar with the elapsed and
  estimated remaining time. Standard output still lists only the files
  written, and the Cube files are identical with and without `--verbose`.

Library changes:

- `ProgressCallback` (`std::function<void(std::size_t done, std::size_t
  total)>`, `<chem/kumo/density.hpp>`) is accepted as an optional last
  parameter by `Molecule::densities_blocked` and `Molecule::orbitals` and as
  `EvaluationOptions::progress`. It is called after each block of grid
  points; calls are serialized (never concurrent), `done` increases up to
  `total`, and an exception thrown by the callback aborts the evaluation
  and reaches the caller. It does not change the results.
- OpenMP: the Evaluator's point path (all density derivatives, and AO-MO
  densities) evaluated its blocks one at a time, so it ran on a single
  thread in OpenMP builds. Its blocks now run on worker threads; each block
  is still passed to the same serial kernel on its own, so results are
  bit-identical to the serial build for any thread count.

Timings against 1.6.0 (HBC def2-TZVP, 1410 AOs, release OpenMP builds,
MacPorts Clang 22, Apple M with 4 performance and 4 efficiency cores; the
absolute numbers depend on the machine):

| Evaluator, 50752 points | 1.6.0, 8 threads | 1.7.0, 8 threads | Speed-up |
|---|---:|---:|---:|
| AO_P gradient | 2.97 s | 0.68 s | ×4.4 |
| AO_P Hessian + Laplacian + tau | 4.06 s | 0.99 s | ×4.1 |
| AO_MO density | 0.43 s | 0.09 s | ×4.8 |
| AO_MO gradient | 1.29 s | 0.30 s | ×4.3 |
| AO_MO Hessian + Laplacian + tau | 3.13 s | 0.73 s | ×4.3 |

With one thread the two versions take the same time.

| `kumo` default grid, 1543260 points | 1.6.0 | 1.7.0 |
|---|---:|---:|
| total density, 1 thread | 17.2 s | 11.9 s |
| total density, 8 threads | 4.9 s | 3.1 s |
| `--property laplacian`, 1 thread | 85.3 s | 76.3 s |
| `--property laplacian`, 8 threads | 85.3 s | 18.9 s |

Validation:

- 71/71 Meson tests passed on macOS ARM64 on 2026-11-12 (serial and
  MacPorts Clang 22 OpenMP builds); CI builds and tests on ubuntu-latest
  and macos-latest. New tests: `progress` (monotonic progress reports,
  documented totals, bit-identical results, exception propagation) and
  `evaluator-parallel` (Evaluator results bit-identical across thread
  counts for AO_P and AO_MO, error propagation from worker blocks). The CLI
  end-to-end test covers `--backend` and `--verbose`. No reference
  tolerance was changed.

---

# chem-kumo 1.6.0 — release notes

Density derivatives from the command line, and faster AO evaluation through
a rigorous per-shell cutoff. The public API is extended in a
source-compatible way; rebuild code that links against libkumo (the layout
of `BasisEvaluator` changed again).

New in the CLI:

- `kumo --property density,gradient,laplacian,tau` writes, for every density
  component chosen with `--spin`, the density (default), the gradient norm
  |grad rho| (e/bohr^4), the Laplacian (e/bohr^5) and the positive kinetic
  energy density tau (hartree/bohr^3) as Cube files. Density files keep
  their names; derivative files are `PREFIX_<component>_<property>.cube`
  (for example `PREFIX_density_laplacian.cube`, `PREFIX_spin_tau.cube`).
  The derivatives use the AO density matrix contraction (AO_P), like the
  density files. `--property` is rejected with `--mo` and `--list-mo`.
  On O2, the integral of the tau Cube matches the kinetic energy derived
  from the FCHK total energy and virial ratio within 0.13% (uniform grid).

Performance:

- Shell cutoff: for each shell, the AO evaluator precomputes a radius
  beyond which a rigorous bound on every AO of the shell and on its first
  and second derivatives is at most `BasisEvaluator::ao_cutoff` = 1e-14; at
  farther points the shell is not computed and gives exact zeros.
  `evaluate_into_reference()` keeps the uncut evaluation. Results change by
  at most about 6e-14 (about 1e-18 relative to the largest value) on O2,
  acrolein and HBC.

| HBC def2-TZVP (1410 AOs), 50752 points, one thread | 1.5.0 | 1.6.0 | Speed-up |
|---|---:|---:|---:|
| `Molecule::density` | 1.05 s | 0.81 s | ×1.3 |
| Evaluator AO_MO density | 0.67 s | 0.43 s | ×1.56 |
| Evaluator AO_P gradient | 3.38 s | 2.98 s | ×1.13 |
| Evaluator AO_P Hessian + Laplacian + tau | 4.50 s | 4.05 s | ×1.11 |
| Evaluator AO_MO gradient | 1.70 s | 1.30 s | ×1.3 |
| Evaluator AO_MO Hessian + Laplacian + tau | 3.60 s | 3.12 s | ×1.15 |

API additions:

- `BasisEvaluator::ao_cutoff` (1e-14), the bound used by the shell cutoff.

Validation:

- 69/69 Meson tests passed on macOS ARM64 on 2026-11-12 (serial and
  MacPorts Clang 22 OpenMP builds); CI builds and tests on ubuntu-latest
  and macos-latest. New test: `ao-shell-cutoff` (every AO, gradient and
  second derivative zeroed by the cutoff is at most `ao_cutoff`, checked
  against the reference evaluator and finite differences). The CLI
  end-to-end test covers `--property`. No reference tolerance was changed.

---

# chem-kumo 1.5.0 — release notes

A performance release: grid densities and density derivatives are several
times faster, most of all for large basis sets. The public API is extended
in a backward-compatible way at source level; rebuild code that links
against libkumo (the layout of `BasisEvaluator` changed).

Timings against 1.4.1 (same program, one thread, release builds, best of
repeated runs, macOS ARM64; the absolute numbers depend on the machine):

| HBC def2-TZVP (1410 AOs), 50752 points | 1.4.1 | 1.5.0 | Speed-up |
|---|---:|---:|---:|
| `Molecule::density` | 10.7 s | 1.10 s | ×9.7 |
| `Molecule::spin_densities` | 10.7 s | 1.69 s | ×6.3 |
| Evaluator AO_P gradient | 22.6 s | 3.54 s | ×6.4 |
| Evaluator AO_P Hessian + Laplacian + tau | 23.9 s | 4.76 s | ×5.0 |
| Evaluator AO_MO density | 1.52 s | 0.69 s | ×2.2 |
| Evaluator AO_MO gradient | 6.55 s | 1.76 s | ×3.7 |
| Evaluator AO_MO Hessian + Laplacian + tau | 13.9 s | 3.81 s | ×3.6 |

On a small molecule (acrolein, 303072 points) the same operations are 1.2
to 2.1 times faster.

Performance changes:

- AO screening (new, on by default): in each block of grid points, AOs
  whose value (and, for derivatives, all their evaluated derivatives) is at
  most `default_screening_threshold` = 1e-12 at every point are left out of
  the contraction, for densities and derivatives with both the AO_P and
  AO_MO backends. The change in the results is at rounding level: at most
  about 1e-15 relative to the largest value in the tested systems (about
  1e-13 e/bohr^3 in the density). A threshold of 0 gives the exact
  contraction.
- Total or spin density alone uses a single density-matrix product; alpha
  and beta densities are contracted directly, and total and difference are
  derived from them, so no request needs more than two products.
- AO derivatives are evaluated in batches with shells in the outer loop, so
  each pure shell needs one Cartesian-to-pure transform per batch and no
  buffers are allocated per point. Results are bit-identical to 1.4.1.
- AO_MO: closed-shell molecules reuse the alpha projection for beta, and
  AO values for density-only requests are evaluated in batches
  (bit-identical results).

API additions:

- `default_screening_threshold` (`<chem/kumo/density.hpp>`) and
  `EvaluationOptions::screening_threshold`; a negative or nonfinite value
  is rejected with `std::invalid_argument`.
- A defaulted `screening_threshold` parameter for
  `Molecule::densities_blocked` and `Molecule::density_derivatives_blocked`.
- `BasisEvaluator::evaluate_derivatives_batch(points, order, output)` and
  `BasisEvaluator::derivative_panel_count(order)` evaluate AO values and
  derivatives for a batch of points into consecutive panels.

Validation:

- 68/68 Meson tests passed on macOS ARM64 on 2026-11-12 (serial and
  MacPorts Clang 22 OpenMP builds); CI builds and tests on ubuntu-latest
  and macos-latest (Homebrew now ships Eigen 5, which is supported). New
  tests: `density-screening` and `ao-derivatives-batch`. No reference
  tolerance was changed.

---

# chem-kumo 1.4.1 — release notes

A maintenance release: stricter and more consistent input validation,
clearer error messages, and internal cleanup. The public API is unchanged
and numerical results are bit-identical to 1.4.0 (checked on O2 and acrolein
outputs, serial and OpenMP builds).

Behaviour changes:

- Grid validation uses one scale-independent rule at every evaluation entry
  point (`Molecule::density*`, `Molecule::orbital*`, `Evaluator::evaluate`):
  positive point counts, finite origin and steps, nonzero steps and linearly
  independent axes, measured on normalized step directions.
  - `Evaluator::evaluate` now rejects degenerate grids (zero steps,
    collinear axes) with `std::invalid_argument`; before, derivative
    requests on such grids returned meaningless values.
  - Orbital evaluation now accepts valid grids with very small spacing,
    which it previously rejected.
  - Error messages for invalid grids are uniform (for example
    "invalid grid geometry: linearly dependent axes").
- `MolecularOrbitals` applies the same checks as `Molecule`: it rejects MO
  records with more MOs than basis functions or nonfinite coefficients.
- A damaged or truncated gzip Cube file is now reported as a gzip error
  (for example `truncated gzip Cube: <path>`) instead of a misleading Cube
  format error. FCHK and Cube readers share one streaming gzip decoder with
  strict CRC32/ISIZE checks; gzip error messages always include the path.

Internal changes:

- Shared internal helpers replace duplicated code for grid points, density
  matrix selection, AO Hessian layout, MO coefficient reading, Cube output
  and OpenMP block loops (first worker exception rethrown on the caller).
- `Molecule::from_fchk` and `Evaluator::evaluate` are split into small,
  documented helpers.
- Code formatted with `.clang-format` and `meson format` (the formatting
  commit is listed in `.git-blame-ignore-revs`); public headers use `///`
  documentation comments.
- Tests report what they expected on failure (`tests/test_support.hpp`).
- CI installs dependencies with retries and timeouts, and installs zlib
  explicitly on Ubuntu.
- 66/66 Meson tests passed on macOS ARM64 on 2026-11-11 (serial and
  MacPorts Clang 22 OpenMP builds); CI builds and tests on ubuntu-latest and
  macos-latest. New tests: `grid-validation` and `parallel-blocks`.

---

# chem-kumo 1.4.0 — release notes

- New typed `chem::kumo::Evaluator` API (`<chem/kumo/evaluator.hpp>`) for
  density components (total, alpha, beta, spin), gradients, Hessians,
  Laplacians, `tau` and orbital amplitudes, driven by an `EvaluationPlan`
  that fuses requests and evaluates only the AO derivatives they need.
- Electronic backends: `Backend::AO_P` (AO density matrix) and
  `Backend::AO_MO` (occupied-orbital contractions). The default
  `Backend::AUTO` selects AO-MO only when the orbital representation stored
  in the FCHK passes validation at load time; an explicitly requested
  ineligible AO-MO backend is rejected. `Results::backend_used()` reports
  the choice. The AO-MO speedup is workload-dependent
  (see `benchmarks/BACKEND_COMPARISON.md`).
- Additional numerical APIs: analytic AO overlap integrals
  (`BasisEvaluator::overlap_matrix()`), Cartesian AO moments through second
  order (`BasisEvaluator::cartesian_moments()`), spin-resolved one-particle
  density matrices (`Molecule::density_matrix()`), metric subspace
  projectors (`SubspaceProjector`), quadrature-based AO weighted overlaps
  (`WeightedOverlapAccumulator`, `weighted_overlap_matrix()`), trilinear
  field sampling (`FieldSampler`) and scalar-field integration
  (`integrate_rectangular()`, `integrate_trapezoidal()`).
- Scalar Gaussian Cube input: `read_cube()` reads plain or gzip-compressed
  Cube files (detected from the gzip signature), decompresses incrementally
  and rejects invalid or truncated streams. Only single-dataset scalar Cube
  files are supported; Cube output remains plain text.
- Atoms now carry an explicit nuclear charge; `Field` samples have
  controlled mutable access; grid geometry is validated independently of
  physical scale.
- All public headers, including the new ones, are installed; the CLI and
  installed-consumer tests derive the expected version from Meson.
- The `kumo` CLI is unchanged: the new facilities are C++ API only.
- 64/64 Meson tests passed on macOS ARM64 on 2026-10-30; CI builds and
  tests on ubuntu-latest and macos-latest. See [docs/API.md](docs/API.md)
  and the README for contracts and validated scope.

---

# chem-kumo 1.3.0 — release notes

- Optional OpenMP grid-block parallelism for density components, spatial
  density derivatives, MO amplitudes and blocked MO spatial derivatives.
- New `orbital_derivatives_blocked()` public API for arbitrary 3 x N points.
- Independent per-block temporary buffers; shared read-only basis and
  coefficient data; errors captured inside workers and rethrown on caller.
- Serial compilation remains the default, without an OpenMP runtime dependency.
- 34/34 functional tests passed on macOS ARM64 (Apple Clang serial and
  MacPorts Clang 22 OpenMP) and Debian x86-64 (GCC 14 OpenMP).
- Benchmark measurements and limitations: `benchmarks/PARALLEL_DENSITY.md`,
  `benchmarks/PARALLEL_FIELDS.md`, and `docs/OPENMP.md`.

---

# chem-kumo 1.2.0 — release notes

This development branch adds analytic Cartesian spatial derivatives of
Gaussian AO basis functions, signed molecular orbitals, and spin-resolved
electron densities. The density derivative API includes the full Hessian,
Laplacian, and positive kinetic-energy density `tau`. Derivatives are
computed from existing FCHK data, with no SCF or nuclear response code.

The density derivative evaluator supports pointwise and Eigen-blocked
paths. The current blocked implementation packs AO values and three
first-derivative panels into one matrix and uses a single dense Eigen
product. Second-derivative panels enter only the linear reductions.
This improves data reuse without changing dense quadratic AO scaling.

An optional spatial benchmark reports AO evaluation, derivative,
pointwise/blocked, and separate/fused contraction timings. On HBC
(B3LYP/def2-TZVP, 1410 AO, 1000 points, Apple Silicon), the measured
blocked B=128 time was 0.4423 s versus 1.3674 s pointwise; the fused
contraction microbenchmark improved by about 4.2%. These are platform-
and workload-specific observations, not guaranteed speedups.

The FCHK reader also accepts files without `P(S=P) Contraction coefficients`
when no SP shell is present. The existing density/MO behavior remains
unchanged. The 31-test suite passed on macOS ARM64 on 2026-10-13, including gzip integrity checks and compressed HBC references.
The CLI does not yet export derivative fields; all new derivative methods
are exposed through the C++ API. Dependencies are C++20, Eigen, and zlib.

See [docs/API.md](docs/API.md),
[mathematical foundations](docs/mathematics/chem-kumo-mathematics.tex), and
[benchmark methodology](benchmarks/SPATIAL_DERIVATIVES.md).

The five FCHK reference datasets are stored as losslessly compressed `.fchk.gz` files, including three HBC benchmarks (678, 720 and 1410 AO). The FCHK reader uses zlib and checks gzip integrity.

---

# chem-kumo 1.1.0 — release notes

This release adds the installable `kumo` command-line tool and molecular-orbital
reconstruction from Gaussian FCHK coefficients. The CLI supports 1-based MO
indices, inclusive ranges, HOMO/LUMO selectors, and alpha/beta channels.
Existing ground-state SCF density functionality is unchanged.

The 25-test development suite passed on macOS ARM64, including ten Gaussian
Cubegen orbital-amplitude comparisons and restricted/unrestricted CLI listing
checks. This release adds an end-to-end CLI regression test for density and
orbital Cube generation and extends the installed-consumer check to cover the
installed CLI. These additions must be validated on the final release package.

The numerical core is unchanged from the validated 1.1.0 development baseline.
Dependencies remain C++20 and Eigen, with Meson as the build system.

OpenMP stage 02: optionally parallelize blocked density grids with a static
schedule, thread-local AO/Eigen matrices, and unchanged arithmetic per block.
Serial builds remain the default; spatial derivatives are still serial.

### OpenMP Stage 03a (development)

- Parallelize `density_derivatives_blocked()` over independent point blocks.
- Preserve fused four-panel density-matrix contractions and point ordering.
- Add serial/OpenMP comparisons for density, gradient, Hessian, Laplacian and tau.
- No changes to orbital evaluation or public API.

### Stage 03c: performance instrumentation

- Add optional `kumo-parallel-fields-benchmark` for density derivatives, one/multiple molecular orbitals, and orbital derivatives.
- Report single-thread versus OpenMP timing, speedup, efficiency, and full-component numerical consistency on matching grids.
- See `benchmarks/PARALLEL_FIELDS.md` for methodology and usage. No numerical kernels are changed.
