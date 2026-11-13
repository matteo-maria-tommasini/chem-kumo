# Optional OpenMP grid evaluation

## Scope and design

OpenMP is opt-in (`-Dopenmp=true`; default `false`). Independent spatial
blocks are assigned with `schedule(static)` in the blocked density,
density-derivative, orbital amplitude and blocked orbital-derivative paths.
Input basis data, density matrices and MO coefficients are shared read-only;
per-block Eigen temporaries are private. Each block writes disjoint output
indices. Exceptions are captured inside worker loops and rethrown after the
parallel region. `EIGEN_DONT_PARALLELIZE` prevents nested Eigen threading.
The `Evaluator` point path (all density derivatives, and AO_MO densities)
runs its blocks of `block_size` grid points on worker threads with
`schedule(dynamic)`; each block is passed to the serial block kernel on its
own, so block boundaries and results do not depend on the thread count.
The derived properties (ELF, reduced gradient, sign(lambda2) rho; 1.8.0)
are computed in the same point path from each block's derivatives. The
electrostatic potential (1.9.0) has its own pass over clusters of
8 x 8 x 8 grid points, also scheduled with `schedule(dynamic)`; cluster
boundaries depend only on the grid and each value only on its point, its
cluster and the far-field plan, so the results do not depend on the thread
count either. The electrostatic potential on a molecular surface (1.10.0)
evaluates the surface vertices in clusters (those on the edges of a block of
8 x 8 x 8 grid cells) in the same way; the mesh itself is built serially,
and its vertices are placed on the isosurface with the blocked density
path, so the results do not depend on the thread count. The pointwise
derivative APIs remain serial. AO screening
(1.5.0) is decided per block of grid points and the AO shell cutoff (1.6.0)
per point, so neither depends on the thread count; see [API.md](API.md).

## Building and testing

```sh
./build-and-test.sh --clean
./build-and-test.sh --clean --openmp --release
```

The helper uses separate `build/` and `build-openmp/` directories. On macOS
it selects the MacPorts Clang 22 native file for OpenMP builds; on Linux
it uses the default toolchain (GCC 14 was tested). A suitable OpenMP runtime
is required only for the OpenMP build. The `libkumo` shared library may depend
on `libomp` (MacPorts) or `libgomp` (GCC); consumers and binary packages must
provide the matching runtime. Check the built library with `otool -L` on
macOS or `ldd` on Linux, and verify installed consumers separately.

Control thread count with `OMP_NUM_THREADS=N`. `OMP_PROC_BIND` and
`OMP_PLACES` can be used for affinity experiments, but no affinity policy
is hard-coded. On dual-socket NUMA machines the optimal policy can vary.

## Validation and benchmarks

The full test suite passes on macOS ARM64 with OpenMP off and on (MacPorts
Clang 22); earlier releases were also checked on Debian x86-64 with GCC 14
and OpenMP enabled. CI (ubuntu-latest, macos-latest, and a sanitizer build
on Ubuntu) runs the serial build only, so OpenMP builds must be tested
locally with `./build-and-test.sh --clean --openmp`. The `parallel-*` tests
compare blocked results across worker counts, and `evaluator-parallel`
requires `Evaluator` results to be bit-identical to one thread for both
backends and several block sizes. Tests and HBC benchmarks found no
difference from one-thread references in the measured cases; this does not
prove universal thread safety.

The optional benchmark targets (`-Dbenchmarks=true`) measure compute paths
without FCHK reading or Cube writing. Times are minimum wall times after a
warm-up; compare speedups only within the same platform and workload.
See [density benchmarks](../benchmarks/PARALLEL_DENSITY.md) and
[field benchmarks](../benchmarks/PARALLEL_FIELDS.md).

The Apple M3 (4P+4E) achieved up to 5.05x on density; dual Xeon E5640
(2 sockets, 8 physical cores / 16 SMT threads) achieved up to 6.89x on
density and 6.77x on density derivatives.

`Evaluator` point path, 1.10.0 release build (MacPorts Clang 22, `-O3`) on
the Apple M3 (4P+4E), HBC with 97,280 grid points (margin 3 bohr, spacing
0.4 bohr), best of three runs, seconds:

| Basis | Request | 1 thread | 4 threads | 8 threads | Speedup (4) | Speedup (8) |
|---|---|---:|---:|---:|---:|---:|
| def2-SVP | gradient | 1.164 | 0.320 | 0.250 | 3.6x | 4.7x |
| def2-SVP | Hessian | 2.688 | 0.819 | 0.632 | 3.3x | 4.3x |
| def2-SVP | AO_MO density | 0.416 | 0.111 | 0.088 | 3.7x | 4.7x |
| def2-TZVP | gradient | 2.430 | 0.724 | 0.575 | 3.4x | 4.2x |
| def2-TZVP | Hessian | 5.779 | 1.770 | 1.386 | 3.3x | 4.2x |
| def2-TZVP | AO_MO density | 0.837 | 0.232 | 0.188 | 3.6x | 4.5x |

Going from 4 to 8 threads adds the four efficiency cores, so the gain is
smaller than from 1 to 4.

These are observed speedups, not performance guarantees. In particular, an OpenMP-enabled library does not
imply that every CLI operation or every pointwise API is parallel.
