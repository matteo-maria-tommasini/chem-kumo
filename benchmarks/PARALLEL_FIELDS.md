# Parallel spatial-derivative and orbital benchmarks (Stage 03c)

This optional benchmark exercises the **public API** without changing numerical
kernels. It measures four independent workloads:

- `density-derivatives`: total-density value, gradient, full Hessian, Laplacian and tau.
- `mo-single`: one occupied alpha MO, evaluated using `Molecule::orbitals()`.
- `mo-multiple`: eight alpha MOs (or fewer if unavailable), evaluated together.
- `mo-derivatives`: one occupied alpha MO value, gradient and full Hessian via
  `Molecule::orbital_derivatives_blocked()`.

All runs use identical coordinates for a given invocation. The density-derivative
and MO-derivative APIs accept arbitrary point matrices; MO amplitude evaluation
uses the matching regular grid. `block_size` affects derivative paths; the
current `orbitals()` implementation has an internal fixed block size of 128.

Parsing the compressed FCHK, constructing the molecule and grid, and generating
point coordinates are excluded from timed sections. Timings **include** field
allocations, evaluation and result collection, but exclude Cube output. Each
thread count receives one untimed warm-up followed by the requested number of
repetitions. The minimum wall time is reported. All scalar output components
are compared to an independently evaluated single-thread reference with an
absolute tolerance of 1e-10; non-finite values are rejected. A zero max_abs
means exact equality in the comparisons performed, not proof of general
thread safety.

Build on macOS with MacPorts Clang 22:

```sh
./build-and-test.sh --clean --openmp --release
meson configure build-openmp -Dbenchmarks=true
meson compile -C build-openmp
```

Example HBC/def2-TZVP, about 10,648 grid points:

```sh
for mode in density-derivatives mo-single mo-multiple mo-derivatives; do
  ./build-openmp/kumo-parallel-fields-benchmark \
    tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz "$mode" 10000 3 128
done
```

Try block sizes 64, 128 and 256 for derivative workloads. The `mo-single`
and `mo-multiple` paths intentionally use their existing internal block size.
A larger grid may be useful for stronger scaling but derivative buffers and
result vectors require substantial RAM. Do not run several large cases
simultaneously. For a controlled comparison, benchmark on an otherwise idle
machine and use the same compiler, build type, grid and MO selections.

OpenMP thread counts are set inside the benchmark to 1, 2, 4, 8 and 16 with
dynamic adjustment disabled. A build without OpenMP runs only the one-thread
reference. Timings should not be compared across machines or toolchains without
accounting for hardware and compiler differences.

## Measured baseline: Apple M3, 2026-10-15

Hardware: Apple M3, 8 physical cores (4 performance + 4 efficiency).
Toolchain: MacPorts Clang 22.1.8, OpenMP 5.1, Meson `release` build.
Input: HBC B3LYP/def2-TZVP (`1410` AO, `60` atoms), `10648` grid points,
`3` timed repetitions per thread count, requested block size `128`.
Times below are **minimum measured wall times** in seconds, not averages.
The MO amplitude paths use their fixed internal block size of 128.

| Workload | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) | Best speedup |
|---|---:|---:|---:|---:|---:|---:|
| Density derivatives | 4.667426 | 2.491336 | 1.364324 | 0.972749 | 0.969306 | 4.815x (16) |
| One MO amplitude | 0.107769 | 0.059484 | 0.031014 | 0.027677 | 0.024157 | 4.461x (16) |
| Eight MO amplitudes | 0.115582 | 0.060875 | 0.033568 | 0.029444 | 0.026856 | 4.304x (16) |
| One MO and spatial derivatives | 0.691053 | 0.418626 | 0.248198 | 0.180457 | 0.174950 | 3.950x (16) |

All reported comparisons gave `max_abs=0` for the components tested.
The machine has only eight physical cores: the 16-thread case is an
oversubscribed experiment, not a recommended default. These are individual
benchmark runs; differences of a few percent should be checked with repeated
independent runs before drawing performance conclusions. The MO derivative
path currently reuses the pointwise analytic evaluator inside its blocked
parallel interface. A possible *future* optimization is batched AO-derivative
contractions; this is **not** implemented in Stage 03.

The validated functional suite passed **34/34** tests with Apple Clang 21
(OpenMP disabled, `debugoptimized`) and **34/34** with MacPorts Clang 22
(OpenMP enabled, `release`). These different build types should not be used
for a direct compiler performance comparison.

## Measured baseline: Avatar, Debian x86-64, 2026-10-15

Dual Intel Xeon E5640, GCC 14.2, OpenMP 4.5, release build. HBC
B3LYP/def2-TZVP, 1410 AO, 5832 points, requested block size 128,
3 timed repetitions per thread count. Minimum wall times in seconds.

| Workload | 1 thread | 2 threads | 4 threads | 8 threads | 16 threads | Best speedup |
|---|---:|---:|---:|---:|---:|---:|
| Density derivatives | 27.959591 | 15.180608 | 9.166416 | 5.850316 | 4.132953 | 6.765x |
| One MO amplitude | 0.514391 | 0.259687 | 0.192623 | 0.139769 | 0.105682 | 4.867x |
| Eight MO amplitudes | 0.547609 | 0.276774 | 0.242195 | 0.122817 | 0.092186 | 5.940x |
| One MO and spatial derivatives | 1.890029 | 0.954644 | 0.744899 | 0.494681 | 0.352898 | 5.356x |

All measured configurations reported `max_abs=0`. MO amplitude timings are
short and susceptible to run-to-run variability; do not interpret small
relative differences as a reproducible optimization.
