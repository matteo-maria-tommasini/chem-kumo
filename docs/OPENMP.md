# Optional OpenMP grid evaluation

## Scope and design

OpenMP is opt-in (`-Dopenmp=true`; default `false`). Without it the same
code runs on one thread.

**How the work is split.** Every parallel loop divides the points into
independent pieces (blocks of consecutive points, or clusters of nearby
points) and gives the pieces to the threads with `schedule(dynamic)`.
Dynamic scheduling matters because the pieces have very different costs:
with AO screening, a block in the empty region far from the molecule is
much cheaper than one near the nuclei, and a static assignment of
contiguous ranges would leave threads idle.

**Why the results do not depend on the threads.** The pieces are fixed by
the input alone (the block size, the grid or the point list), never by the
number of threads; each piece is computed on its own and writes its own
output entries. The arithmetic of every value is therefore the same for any
number of threads, and results are identical bit for bit. Basis data,
density matrices and MO coefficients are shared read-only; each thread owns
its temporary matrices. An exception thrown on a thread is caught there and
rethrown on the calling thread after the loop. `EIGEN_DONT_PARALLELIZE`
prevents Eigen from starting threads of its own inside these loops.

**Parallel loops.**

| Computation | Pieces |
|---|---|
| Densities on grids (`densities_blocked`), orbital amplitudes on grids | blocks of grid points |
| Density derivatives and AO-MO densities (`Evaluator` point path), with ELF, LOL, the reduced gradient and the other derived fields | blocks of `block_size` points |
| Blocked density and orbital derivatives (`density_derivatives_blocked`, `orbital_derivatives_blocked`) | blocks of points |
| Orbital-energy fields (ALIE, LEA, LEAE, local electronegativity and hardness) and orbitals at arbitrary points (1.13.0) | blocks of `block_size` points |
| Electrostatic potential on grids (1.9.0) | tiles of 8 x 8 x 8 grid points |
| Electrostatic potential on isosurface vertices (1.10.0) | vertices on the edges of a block of 8 x 8 x 8 grid cells |
| Electrostatic potential at arbitrary points (1.13.0) | points in the same cube of 1.6 bohr |
| ALIE on a surface (`local_ionization_energy_surface`) | blocks of 128 vertices |

The isosurface mesh is built serially, and its vertices are placed on the
surface with the blocked density path. Mapping fields on a surface
(`map_on_density_isosurface`, 1.13.0) evaluates the vertices with the
`Evaluator` at points, so it uses the loops above. The pointwise derivative
functions (`density_derivatives`, `orbital_derivatives` at one point) are
serial. AO screening is decided per block of points and the AO shell
cutoff per point, so neither depends on the number of threads; see
[API.md](API.md).

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

Blocked density (`densities_blocked`), release build on the Apple M3
(4P+4E), HBC with 103,823 points and block size 256, with
`schedule(dynamic)`: 3.7x (def2-SVP) and 3.5x (def2-TZVP) with 4 threads,
5.0x and 4.8x with 8 threads. The 1.10.0 release, with `schedule(static)`,
gave only 2.8x and 2.5x with 4 threads (4.2x and 3.8x with 8), because
screened blocks far from the molecule are much cheaper and contiguous
block ranges left threads idle; see
[density benchmarks](../benchmarks/PARALLEL_DENSITY.md). Before AO screening
(1.5.0) the same benchmark gave up to 5.05x on the M3, and a dual Xeon E5640
(2 sockets, 8 physical cores / 16 SMT threads) up to 6.89x on density and
6.77x on density derivatives; the Xeon has not been measured again.

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
