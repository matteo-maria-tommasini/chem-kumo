# AO–P versus occupied-MO benchmark (Stage 01A/01B)

This is a **single-thread, density-only, diagnostic benchmark**. It does not change the chem-kumo API or production kernels.

Build with `meson configure build-openmp -Dbenchmarks=true` and `meson compile -C build-openmp`.

Run:

```sh
OMP_NUM_THREADS=1 ./build-openmp/kumo-backend-comparison tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz 4096 5 128
OMP_NUM_THREADS=1 ./build-openmp/kumo-backend-comparison tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz 4096 5 128
```

The tool reads `Total SCF Density` from the FCHK, unpacks its triangular AO matrix, and separately reads occupied alpha/beta MO coefficients. For a restricted wavefunction without beta coefficients, it uses twice the alpha orbital density. For unrestricted wavefunctions it sums occupied alpha and beta orbital densities. This assumes **unit occupations for the first `Nalpha`/`Nbeta` orbitals** and does not yet handle fractional occupations. Differences from the FCHK density are *reported*, not suppressed. Large discrepancies indicate that the orbital representation may not reproduce the stored density and invalidate direct scientific equivalence.

Measurements: `ao_time_s` evaluates and stores AO values on the full grid; `density_matrix_time_s` and `occupied_orbitals_time_s` time contractions with **already evaluated AO values**. Each measurement is warmed up, and the minimum of the requested repetitions is reported. Parsing and matrix preparation are excluded. FLOP counts are nominal leading-order dense GEMM estimates, not hardware counters. The density-only AO–P path uses one P*Phi product; occupied-MO uses one or two C^T*Phi products.

**Limitations:** This stage does not yet benchmark derivatives, true fused end-to-end execution, memory allocations, or OpenMP scaling. The AO timing currently includes copying into a full AO cache, which is needed for a fair shared-input kernel comparison but is not representative of streaming production execution. The full AO cache costs `8*N*G` bytes; keep point counts moderate. The two methods can differ due to floating-point rounding or FCHK data conventions. Reported relative error uses a denominator floor of `1e-10`.

Future stages should add a fused streaming end-to-end benchmark, shared preallocated per-worker scratch buffers, derivatives, and a measured backend crossover.

## Stage 01B: density-matrix consistency and streaming

The additional `density_matrix_reconstruction_*` metrics compare the FCHK
`Total SCF Density` matrix directly with the AO density reconstructed from
integer-occupied MO coefficients. The maximum absolute, maximum relative
(normalized by the largest absolute FCHK matrix entry), and relative Frobenius
norm differences are reported. These checks do **not** assume that every FCHK
represents an integer-occupation SCF density; discrepancies must be interpreted
before selecting an MO backend.

The `streaming_*` timings perform AO evaluation and density construction in
one block loop, reusing allocated AO, contraction, and output buffers. They do
not read from the full-grid AO cache used by the kernel-only measurements.
The two streaming paths are timed independently with one thread. Parsing,
coordinate generation, and matrix reconstruction are excluded. The minimum
of the requested repeats is reported, after a warm-up. No screening or
approximations are introduced. `streaming_max_abs` compares the two streamed
densities. This is a diagnostic prototype, **not** the production Evaluator.

Caution: the benchmark still creates the full-grid AO cache for the separate
kernel-only timings, even though the streaming measurements do not access it.
The reported times are not directly additive: kernel-only AO cache timings
include writes that the streaming path avoids. Benchmark runs may vary due to
CPU frequency, memory contention, and compiler/runtime configuration.

## Stage 01D: component-resolved numerical diagnostics

The derivative benchmark additionally prints `component=` records for each
computed scalar, gradient component, kinetic-energy density, and Hessian
component. `max_abs` and `max_rel_floor_1e-10` compare the AO-density-matrix
and occupied-MO implementations; `worst_point` is the zero-based grid index
of the largest absolute discrepancy. The relative metric uses the AO-P value
with a fixed absolute denominator floor of `1e-10`, so it must not be
interpreted as a meaningful relative error in the density tail.

For `density-hessian`, it constructs symmetric 3x3 Hessians from the six
components and compares **ordered** eigenvalues. `hessian_eigen_bound_excess`
checks a conservative Frobenius-norm version of Weyl's eigenvalue bound;
`eigen_trace_max_abs` checks that the sum of eigenvalues equals the matrix
trace, and `laplacian_max_abs` compares traces of the two Hessians. These
trace checks are algebraic consistency diagnostics, **not** an independent
Laplacian calculation. `lambda2_sign_ambiguous` counts points where the
middle eigenvalue of either representation lies within the Frobenius-norm
perturbation bound; these are not reliable sign classifications.
`lambda2_sign_stable_mismatch` must remain zero.

`low_density_threshold` reports how many sampled points have AO-P density
magnitude below each diagnostic cutoff (atomic units). These counts do not
change or truncate any density and are **not** ELF/RDG numerical policies.

This benchmark uses integer occupations and compares representations derived
from finite-precision FCHK data. It is a diagnostic, not a replacement for
independent scientific reference tests. The derivative AO evaluator is still
pointwise, and the benchmark does not measure OpenMP scaling.
