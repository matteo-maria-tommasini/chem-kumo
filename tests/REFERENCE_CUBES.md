# Gaussian Cube reference files

This directory contains Gaussian 09 reference Cube files used for
quantitative validation of chem-kumo.

The electron-density reference files were generated using the
full SCF electron density:

    gcubegen09d 0 FDensity=SCF o2.fchk o2_density.cube -3 h
    gcubegen09d 0 FDensity=SCF acrolein.fchk acrolein_density.cube -3 h

The corresponding numerical electron counts on the reference grids
are approximately:

| System | Integrated electron density |
|--------|-----------------------------|
| O2 | 16.0872 |
| Acrolein | 30.3105 |

These values include both core and valence electrons.

The deviations from the exact electron counts (16 and 30) arise
from numerical integration on finite three-dimensional grids.

The automated reference tests compare chem-kumo results against
Gaussian Cube data point by point, including numerical integration
and error statistics.

The corresponding FCHK inputs are stored as `o2.fchk.gz` and `acrolein.fchk.gz`.

## Multiwfn references for derived properties

`tests/data/multiwfn/` contains Cube files of the electron localization
function (ELF), reduced density gradient (RDG) and sign(lambda2) rho,
computed with Multiwfn 2026.9.1 (`Multiwfn_noGUI`, macOS build) on
2026-11-12. The inputs are the FCHK files of this repository:

| File | FCHK | Multiwfn function |
|------|------|-------------------|
| `h2o_elf.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 9, ELF |
| `h2o_rdg.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 13, RDG |
| `h2o_sign_lambda2_rho.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 15, sign(lambda2)*rho |
| `o2_elf.cube.gz` | `o2.fchk.gz` | 9, ELF (open shell) |

**Settings.** Multiwfn ran with its distributed `settings.ini`, with two
changes so that it computes the plain definitions that chem-kumo
implements:

- `ELF_addminimal= 0`. The default 1 adds 1e-5 to D(r).
- `RDG_maxrho= 0`. The default 0.05 sets RDG to 100 where rho >= 0.05.

The other relevant settings keep their default values: `ELFLOL_type= 0`
(Becke's ELF) and `expcutoff= -40`.

**Grid and output.** The grid spans the molecule with an extension of 3
bohr and a spacing of 0.3 bohr. The output was exported as a Gaussian-type
Cube file and compressed with `gzip -9 -n`. For example, the water ELF was
produced with this input on standard input:

    h2o.fchk
    5
    9
    -10
    3
    4
    0.3
    2
    0
    q

where 5 calculates grid data, 9 selects ELF, `-10` and `3` set the grid
extension to 3 bohr, `4` and `0.3` set the grid spacing, 2 exports the
Cube file, then 0 and q leave the program. The other files use functions
13 and 15, or `o2.fchk`.

**Test.** The test `multiwfn-reference` evaluates the same quantities with
the Evaluator on the grid of each Cube header and compares them point by
point. The Cube files store six significant digits, so rounding alone
allows a relative difference of 5e-6. The test accepts 2e-5 relative. The
largest differences observed are:

| Quantity | Max relative difference |
|----------|-------------------------|
| Water ELF | 5.9e-6 |
| Water RDG | 1.2e-5, near the O-H bond critical point, where the gradient nearly cancels |
| Water sign(lambda2) rho | 9.0e-6, at the oxygen nucleus |
| O2 ELF | 7.0e-6 |

The O2 file confirms that both programs use the spin-polarized form of the
total ELF for open shells. There are no sign differences in sign(lambda2)
rho.

**Citation.** Multiwfn asks that work using it cite T. Lu and F. Chen,
J. Comput. Chem. 33, 580 (2012) and T. Lu, J. Chem. Phys. 161, 082503
(2024).
