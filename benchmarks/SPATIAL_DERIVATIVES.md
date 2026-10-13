# Spatial derivative benchmark (developer tool)

`kumo-spatial-benchmark` measures the cost of evaluating spatial derivatives
from a **fixed, converged FCHK**. It does not run Gaussian or perform SCF.
The benchmark is optional and does not run as part of the numerical CI tests.

## Build and run

```sh
./build-and-test.sh --clean
meson configure build -Dbenchmarks=true -Dbuildtype=release
meson compile -C build
./build/kumo-spatial-benchmark tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz 1000 5
```

Usage: `kumo-spatial-benchmark FILE.fchk[.gz] [points=125] [repeats=2] [margin_bohr=2.5]`.
The actual point count is the smallest cubic lattice (at least two points
per side) containing at least the requested count. Coordinates span the
nuclear bounding box plus the margin, in bohr. For cross-basis comparisons,
use the **same geometry, point count, and margin**; verify the printed point
count. The benchmark performs one warm-up and reports the *minimum* elapsed
wall time across repeats. Use an otherwise idle machine and repeat complete
runs to estimate variability.

## Output fields and interpretation

| Field | Meaning |
| --- | --- |
| `basis`, `atoms`, `points`, `repeats` | Problem dimensions and sampling |
| `ao_values_s` | AO value evaluation over the selected points |
| `ao_derivatives_s` | AO values, gradients, and full Hessians |
| `ao_derivatives_over_values` | Relative cost of AO derivatives |
| `contraction_block_points` | Number of columns in the separate contraction microbenchmark (up to 64) |
| `four_dense_contractions_s` | Four separate `P * panel` products |
| `fused_dense_contraction_s` | One `P * [panel0 panelx panely panelz]` product |
| `fused_contraction_speedup` | Separate/fused time ratio for the microbenchmark |
| `fused_contraction_max_abs` | Maximum absolute difference between the two product layouts |
| `density_only_s` | End-to-end ordinary density evaluation |
| `derivatives_pointwise_s` | End-to-end pointwise density, gradient, Hessian, Laplacian and tau |
| `derivatives_block_{16,32,64,128,256}_s` | End-to-end blocked derivatives for each block size |
| `speedup` | Pointwise time divided by corresponding blocked time |
| `checksum` | A lightweight guard against discarding computed results, **not** a validation metric |

The four separate products and the fused product use the same AO panels and
selected density matrix. Their times exclude AO evaluation and most of the
end-to-end allocation/reduction work. They are **not** additive components
of the end-to-end measurements. `fused_contraction_max_abs=0` is an observed
numerical comparison for the tested input, not a mathematical guarantee of
bitwise equality on every architecture.

The blocked path computes AO derivatives point by point, packs values and
first derivatives into a contiguous Eigen panel, performs one dense product
with the AO density matrix, and reduces the six Hessian panels without
additional matrix products. This changes data organization, **not** the
$O(N_AO^2 N_points)$ dense contraction scaling. The public API currently
returns one `DensityDerivatives` object per point; it does not stream fields.

## Measured Apple Silicon example (2026-10-13)

HBC (60 atoms), B3LYP/def2-TZVP, **1410 AO**, 1000 grid points, five
repetitions, Apple Clang 21, Eigen 3.4.1, Meson release build:

| Measurement | Before fused product | After fused product |
| --- | ---: | ---: |
| Pointwise derivatives | 1.378945 s | 1.367380 s |
| Blocked derivatives, B=64 | 0.468500 s | 0.447480 s |
| Blocked derivatives, B=128 | 0.460114 s | 0.442308 s |
| Blocked derivatives, B=256 | 0.451660 s | 0.442346 s |

The direct microbenchmark on the fused build reported
`four_dense_contractions_s=0.024762708`,
`fused_dense_contraction_s=0.023753916`, and
`fused_contraction_speedup=1.04246845025` (about **4.2%**).
The blocked B=128 path was about **3.09x** faster than the pointwise path
in that run. Differences between separate complete runs are not controlled
A/B measurements and can include thermal, scheduling, and cache effects.
No speedup is guaranteed for other systems or machines.

Three HBC FCHK references are committed in compressed form under
`tests/data/hbc/`: B3LYP/6-31G(d,p) (720 AO), B3LYP/def2-SVP
(678 AO), and B3LYP/def2-TZVP (1410 AO). All have 60 atoms and 270
electrons. The normal numerical test suite checks their FCHK headers and
compressed integrity; the full HBC timing benchmarks remain opt-in.
There are no timing thresholds in Meson tests.

## Next optimization boundary

Spatial AO screening is **not implemented**. Any future screening must
account for AO values, gradients and Hessians and must be validated against
the unscreened dense path before becoming a default. The present code uses
only Eigen as a numerical dependency.
