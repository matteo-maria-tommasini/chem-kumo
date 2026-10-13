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
