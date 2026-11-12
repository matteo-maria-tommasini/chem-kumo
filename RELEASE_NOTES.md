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
