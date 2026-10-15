# Parallel density benchmark (OpenMP Stage 02b)

This is a **developer benchmark**, not a correctness or timing-gated CI test.
It exercises `Molecule::densities_blocked()` (total and spin-difference
fields) on a fixed grid for a single FCHK, sweeping **1, 2, 4, 8, 16**
OpenMP threads in the same executable. With OpenMP disabled it measures only
one thread. The grid and block size remain identical throughout the sweep.

## Build on macOS

```sh
./build-and-test.sh --clean --openmp --release
meson configure build-openmp -Dbenchmarks=true
meson compile -C build-openmp
```

## Build on Debian / GCC

```sh
./build-and-test.sh --clean --openmp --release
meson configure build-openmp -Dbenchmarks=true
meson compile -C build-openmp
```

## Run

```sh
./build-openmp/kumo-parallel-density-benchmark \
    tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz 10000 3 128
```

Arguments: `FILE [points=10000] [repeats=3] [block_size=128]
[margin_bohr=2.5]`. The actual number of points is `ceil(cbrt(points))^3`
(e.g. 10648 for 10000 requested points). For block-size comparisons,
repeat with 64 and 256; do not change the grid between runs.

Each measurement includes allocation of output fields and all density
contractions, but **excludes** FCHK parsing, molecular construction, grid
construction and Cube writing. One untimed warm-up is performed per thread
count, then the **minimum wall time** of the requested repeats is reported.
The reference is evaluated with one OpenMP thread. Every thread count is
checked against the reference for both output fields (`max_abs <= 1e-10`).
A checksum is printed outside the timed region.

`speedup` is relative to the one-thread timing of **the same binary**;
`efficiency = speedup / threads`. Requesting 16 threads does not imply the
machine has 16 physical cores; oversubscription and heterogeneous cores can
reduce scaling. Run on an otherwise idle machine. For a fair comparison of
OpenMP and non-OpenMP binaries, build both with the **same compiler and
optimization flags**. This benchmark does not measure CLI end-to-end time.

## Measured scaling: Apple M3 (2026-10-15)

Machine: Apple M3, macOS ARM64, 8 physical cores (4 performance + 4
efficiency), MacPorts Clang 22.1.8, OpenMP 5.1. Release build. These are
observed wall times, not a performance guarantee. Each timing is the minimum
of three measured runs after an untimed warm-up for that thread count.
The benchmark evaluates total and spin-difference densities and excludes
FCHK loading, grid construction and Cube output.

### HBC / B3LYP / def2-TZVP, 1410 AOs, 10,648 points

| Block size | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) |
|---:|---:|---:|---:|---:|---:|
| 64 | 2.286467 | 1.153653 | 0.602648 | 0.473768 | 0.476802 |
| 128 | 2.172488 | 1.106114 | 0.578708 | 0.451698 | 0.462691 |
| 256 | 2.099165 | 1.078350 | 0.596394 | 0.431852 | 0.427638 |

### HBC / B3LYP, 103,823 points, block size 256

| Basis | AOs | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) |
|---|---:|---:|---:|---:|---:|---:|
| def2-SVP | 678 | 4.948161 | 2.543645 | 1.332847 | 0.992765 | 0.979464 |
| def2-TZVP | 1410 | 20.478881 | 10.410269 | 5.663391 | 4.244277 | 4.312944 |

Speedup relative to one thread with the **same basis and block size**:

| Basis | 2 threads | 4 threads | 8 threads | 16 threads |
|---|---:|---:|---:|---:|
| def2-SVP | 1.9453 | 3.7125 | 4.9842 | 5.0519 |
| def2-TZVP | 1.9672 | 3.6160 | 4.8251 | 4.7482 |

All reported configurations returned `max_abs=0` relative to the one-thread
reference for both computed fields. The near-ideal scaling up to four threads
and the plateau around eight threads are consistent with, but do not prove,
effects from the M3's heterogeneous cores and/or memory hierarchy. Increasing
the grid from 10,648 to 103,823 points did not improve the eight-thread
speedup for def2-TZVP. Differences of around one percent between 8 and 16
threads should not be interpreted as significant without more independent runs.

Reproduce the larger-grid measurements with:

```sh
./build-openmp/kumo-parallel-density-benchmark \
    tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz 100000 3 256
./build-openmp/kumo-parallel-density-benchmark \
    tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz 100000 3 256
```

These timings measure only the parallelized `densities_blocked()` path.
They are not directly comparable with `gcubegen09d` end-to-end timings or
with a different compiler, grid, or set of requested fields.
