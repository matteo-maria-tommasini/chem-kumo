# Optional OpenMP grid evaluation

## Scope and design

OpenMP is opt-in (`-Dopenmp=true`; default `false`). Independent spatial
blocks are assigned with `schedule(static)` in the blocked density,
density-derivative, orbital amplitude and blocked orbital-derivative paths.
Input basis data, density matrices and MO coefficients are shared read-only;
per-block Eigen temporaries are private. Each block writes disjoint output
indices. Exceptions are captured inside worker loops and rethrown after the
parallel region. `EIGEN_DONT_PARALLELIZE` prevents nested Eigen threading.
The pointwise derivative APIs remain serial. No AO screening is introduced.

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

The 34-test suite passed on macOS ARM64 with OpenMP off/on, and on Debian
x86-64 with GCC 14 OpenMP enabled. Tests and HBC benchmarks compared
parallel results to one-thread references with reported `max_abs=0` in the
measured cases; this does not prove universal thread safety.

The optional benchmark targets (`-Dbenchmarks=true`) measure compute paths
without FCHK reading or Cube writing. Times are minimum wall times after a
warm-up; compare speedups only within the same platform and workload.
See [density benchmarks](../benchmarks/PARALLEL_DENSITY.md) and
[field benchmarks](../benchmarks/PARALLEL_FIELDS.md).

The Apple M3 (4P+4E) achieved up to 5.05x on density; dual Xeon E5640
(2 sockets, 8 physical cores / 16 SMT threads) achieved up to 6.89x on
density and 6.77x on density derivatives. These are observed speedups, not
performance guarantees. In particular, an OpenMP-enabled library does not
imply that every CLI operation or every pointwise API is parallel.
