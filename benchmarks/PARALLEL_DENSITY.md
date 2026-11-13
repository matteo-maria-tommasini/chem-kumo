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

## Measured scaling: Apple M3, dynamic schedule (2026-11-13)

Machine: Apple M3, macOS 27.2 ARM64, 8 physical cores (4 performance + 4
efficiency), MacPorts Clang 22.1.8, release build (`-O3`) of 1.10.0 with the
blocked loops changed from `schedule(static)` to `schedule(dynamic)`, with AO
screening (1.5.0) and the AO shell cutoff (1.6.0). These are observed wall
times, not a performance guarantee. Each timing is the minimum of three
measured runs after an untimed warm-up for that thread count. All reported
configurations returned `max_abs=0` relative to the one-thread reference,
and the checksums equal those of the static schedule.

### HBC / B3LYP / def2-TZVP, 1410 AOs, 10,648 points

| Block size | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.322884 | 0.169887 | 0.096825 | 0.078128 | 0.078711 |
| 128 | 0.385347 | 0.203386 | 0.111042 | 0.086768 | 0.085437 |
| 256 | 0.541672 | 0.284856 | 0.163322 | 0.122702 | 0.120586 |

With screening, smaller blocks are faster: an AO is kept for a block if it
is significant at any point of the block, so smaller blocks keep fewer AOs.

### HBC / B3LYP, 103,823 points, block size 256

| Basis | AOs | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) |
|---|---:|---:|---:|---:|---:|---:|
| def2-SVP | 678 | 0.796900 | 0.404577 | 0.214061 | 0.158176 | 0.152267 |
| def2-TZVP | 1410 | 2.687340 | 1.405087 | 0.759516 | 0.561162 | 0.580688 |

Speedup relative to one thread with the **same basis and block size**:

| Basis | 2 threads | 4 threads | 8 threads | 16 threads |
|---|---:|---:|---:|---:|
| def2-SVP | 1.97 | 3.72 | 5.04 | 5.23 |
| def2-TZVP | 1.91 | 3.54 | 4.79 | 4.63 |

A first A/B run of the same change gave 3.67x / 4.24x (def2-SVP) and
3.59x / 4.78x (def2-TZVP) with 4 / 8 threads, so the four-thread gain is
reproducible, while the eight-thread def2-SVP figure varied between runs.
Differences of a few percent between 8 and 16 threads should not be
interpreted as significant without more independent runs.

## Static schedule: Apple M3, release 1.10.0 (2026-11-13)

Same machine and build as above, but with the blocked loops still using
`schedule(static)`, as released in 1.10.0.

### HBC / B3LYP / def2-TZVP, 1410 AOs, 10,648 points

| Block size | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.322132 | 0.171428 | 0.129271 | 0.092294 | 0.080945 |
| 128 | 0.392738 | 0.206097 | 0.155280 | 0.120261 | 0.095561 |
| 256 | 0.547389 | 0.287471 | 0.237245 | 0.144920 | 0.127461 |

### HBC / B3LYP, 103,823 points, block size 256

| Basis | AOs | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 16 threads (s) |
|---|---:|---:|---:|---:|---:|---:|
| def2-SVP | 678 | 0.802822 | 0.412830 | 0.285763 | 0.189389 | 0.175049 |
| def2-TZVP | 1410 | 2.719077 | 1.438418 | 1.066839 | 0.723465 | 0.619517 |

| Basis | 2 threads | 4 threads | 8 threads | 16 threads |
|---|---:|---:|---:|---:|
| def2-SVP | 1.94 | 2.81 | 4.24 | 4.59 |
| def2-TZVP | 1.89 | 2.55 | 3.76 | 4.39 |

Screening made the one-thread times 6x (def2-SVP) and 7.5x (def2-TZVP)
shorter than in the 2026-10-15 measurements below, but with the static
schedule the four-thread speedup dropped from about 3.7x to 2.5-2.8x; two
further def2-SVP sweeps gave the same 2.8x. The cause was load imbalance,
confirmed by the dynamic-schedule measurements above: `schedule(static)`
gives each thread a contiguous range of blocks, and with screening the
blocks in the empty outer region of the grid cost much less than those
near the molecule. Dynamic scheduling restores the four-thread speedup to
about 3.5-3.7x without changing any result.

## Earlier measurement: Apple M3, before AO screening (2026-10-15)

Same machine, MacPorts Clang 22.1.8, OpenMP 5.1, release build, before AO
screening (1.5.0) and the shell cutoff (1.6.0).

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

| Basis | 2 threads | 4 threads | 8 threads | 16 threads |
|---|---:|---:|---:|---:|
| def2-SVP | 1.9453 | 3.7125 | 4.9842 | 5.0519 |
| def2-TZVP | 1.9672 | 3.6160 | 4.8251 | 4.7482 |

All reported configurations returned `max_abs=0`.

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

## Measured baseline: Avatar, Debian x86-64, 2026-10-15 (before AO screening)

Dual Intel Xeon E5640 (2 sockets, 4 physical cores/socket, 2 SMT threads/core,
2 NUMA nodes), GCC 14.2, OpenMP 4.5, release build. HBC B3LYP/def2-TZVP,
1410 AO, 5832 points, 3 timed repetitions; minimum wall times in seconds.

| Threads | Block 64 | Block 128 | Block 256 |
|---:|---:|---:|---:|
| 1 | 14.312816 | 14.010070 | 13.859962 |
| 2 | 7.280535 | 7.144165 | 7.175005 |
| 4 | 4.524861 | 4.545798 | 4.511535 |
| 8 | 3.011630 | 2.886840 | 2.846716 |
| 16 | 2.207877 | 2.033396 | 2.272181 |

At block 128, the 16-thread speedup was 6.88999x relative to the
same-block one-thread timing; all reported `max_abs=0`. With block 128,
`OMP_PROC_BIND=close OMP_PLACES=cores` gave 2.130629 s at 16 threads,
and `spread` gave 2.223092 s, versus 2.033396 s with default placement.
These are single benchmark sweeps; small differences require independent
replication before drawing conclusions about NUMA affinity.
