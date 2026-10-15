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
