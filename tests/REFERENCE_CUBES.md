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
function (ELF), reduced density gradient (RDG), sign(lambda2) rho and the
electrostatic potential (ESP), computed with Multiwfn 2026.9.1 (`Multiwfn_noGUI`, macOS build) on
2026-11-12. The inputs are the FCHK files of this repository:

| File | FCHK | Multiwfn function |
|------|------|-------------------|
| `h2o_elf.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 9, ELF |
| `h2o_rdg.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 13, RDG |
| `h2o_sign_lambda2_rho.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 15, sign(lambda2)*rho |
| `o2_elf.cube.gz` | `o2.fchk.gz` | 9, ELF (open shell) |
| `h2o_esp.cube.gz` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` | 12, total ESP |
| `o2_esp.cube.gz` | `o2.fchk.gz` | 12, total ESP (open shell) |

The ESP files were added later the same day.

**Settings.** Multiwfn ran with its distributed `settings.ini`, with two
changes so that it computes the plain definitions that chem-kumo
implements:

- `ELF_addminimal= 0`. The default 1 adds 1e-5 to D(r).
- `RDG_maxrho= 0`. The default 0.05 sets RDG to 100 where rho >= 0.05.
- `iESPcode= 1`, for the ESP files: Multiwfn's original analytic ESP code.
  The default 2 (the LIBRETA library, fast version) differed from both the
  original code and chem-kumo by up to 3e-7 hartree/e on the water grid;
  the original code agrees with chem-kumo within the rounding of the Cube
  file.

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
13, 15 and 12 (ESP), or `o2.fchk`.

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

For the ESP the test accepts 2e-5 relative plus 1e-7 hartree/e absolute.
The grid points are read from the Cube header, whose six decimals move them
by up to 5e-7 bohr, and the potential changes by up to about 5e-8 over that
distance away from the nuclei. The largest differences are 0.49 (water) and
0.32 (O2) of the allowed difference; where |V| >= 0.01 the relative
difference is at most 1e-5. The O2 grid contains one nucleus, where
Multiwfn writes about 2e16 and chem-kumo the potential at the nucleus; the
test skips points within 1e-3 bohr of a nucleus.

The O2 file confirms that both programs use the spin-polarized form of the
total ELF for open shells. There are no sign differences in sign(lambda2)
rho.

### Density fields of 1.12.0

Cube files of the fields added in 1.12.0 were computed on 2026-11-13 with
the same Multiwfn build, grids (extension 3 bohr, spacing 0.3 bohr) and
input sequence, from a scratch copy of the same `settings.ini` with these
further settings, so that Multiwfn computes the plain definitions:

- `IRI_rhocut= 0`. The default 5e-5 sets IRI to 5 where rho <= 5e-5.
- `ipolarpara= 1` for the spin polarization: function 5 then returns
  (rho_alpha - rho_beta) / (rho_alpha + rho_beta) instead of the spin
  density.
- `iuserfunc` selects Multiwfn's user-defined function 100.

| File | FCHK | Multiwfn function |
|------|------|-------------------|
| `h2o_lol.cube.gz` | water | 10, LOL |
| `o2_lol.cube.gz` | O2 | 10, LOL (open shell) |
| `h2o_iri.cube.gz` | water | 24, IRI |
| `h2o_dori.cube.gz` | water | 100 with `iuserfunc= 20`, DORI |
| `h2o_kinetic_k.cube.gz` | water | 6, K(r) |
| `o2_kinetic_k.cube.gz` | O2 | 6, K(r) |
| `h2o_potential_energy.cube.gz` | water | 100 with `iuserfunc= 10`, V(r) |
| `h2o_energy_density.cube.gz` | water | 100 with `iuserfunc= 11`, H(r) |
| `o2_spin_polarization.cube.gz` | O2 | 5 with `ipolarpara= 1` |
| `h2o_shannon.cube.gz` | water | 100 with `iuserfunc= 50`, -rho ln rho |
| `h2o_fisher.cube.gz` | water | 100 with `iuserfunc= 51`, \|grad rho\|^2 / rho |

Multiwfn computes K as -(1/2) sum_i psi_i lap psi_i, V as -K - G and H as
-K (`function.f90`), the same quantities as kumo's
`tau - lap rho / 4` forms.

**Test.** `multiwfn-reference` compares these files like the others. Some
fields are small differences of large terms, or small where a ratio
vanishes, so a relative tolerance on the value alone would test rounding
in the cancellation rather than the field. The tolerances, about three
times the largest observed differences, are:

| Field | Tolerance | Largest difference |
|-------|-----------|--------------------|
| LOL (water, O2) | 2e-5 relative | 5.2e-6 relative |
| IRI | 5e-5 relative | 1.5e-5, near the O-H bond critical point (as RDG) |
| DORI | 2e-5 relative + 1e-5 | 6.5e-6 absolute, where theta is about 0.4 |
| K, V, H (water, O2) | 5e-5 of \|tau\| + \|lap rho\| / 4 | 1.6e-5 of that scale |
| Spin polarization (O2) | 2e-5 relative + 1e-6 | 5.8e-7 absolute (Cube rounding at 0.12) |
| Shannon density | 2e-5 of rho (1 + \|ln rho\|) | 5.9e-6 of that scale, near rho = 1 |
| Fisher density | 1e-4 relative | 3.4e-5, near the O-H bond critical point |

The test reports each largest difference as a fraction of its tolerance; it
is at most 0.41.

### Lines and planes (test `line-plane`)

The test `line-plane` reuses `h2o_lol.cube.gz` and `h2o_elf.cube.gz`, with
no new files. `Grid::line` and `Grid::plane` are built from corner points of
the Cube grid, with margin 0 and the Cube spacing (0.3 bohr), so that they
pass through its points: the line is the row j = 10, k = 10 (30 points
along x) and the plane the slice j = 10 (y = 0, the molecular plane; 30 x
24 points). The points agree with the Cube points to 1e-12 bohr, and LOL
and ELF with the Multiwfn values within the 2e-5 relative tolerance used
above; the largest difference is 0.24 of it.

## Multiwfn references for the ESP surface analysis

`tests/data/multiwfn/` also contains the output of Multiwfn's quantitative
molecular surface analysis (main function 12) of the ESP on the
rho = 0.001 e/bohr^3 isosurface of the total density, computed with the
same Multiwfn build and settings (`iESPcode= 1`) on 2026-11-13:

| File | FCHK |
|------|------|
| `h2o_esp_surface.txt` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` |
| `o2_esp_surface.txt` | `o2.fchk.gz` |
| `acrolein_esp_surface.txt` | `acrolein.fchk.gz` |

Each file is the part of the output from the global extrema to the end of
the summary, with terminal colour codes and Multiwfn's floating-point
exception note removed. The water analysis was run with this input:

    h2o.fchk
    12
    3
    0.05
    0
    -1
    -1
    q

where 12 opens the surface analysis (default surface: rho = 0.001; mapped
function: ESP), 3 and 0.05 set the grid spacing to 0.05 bohr instead of
the default 0.25, 0 starts the analysis, and the rest leaves the program.
All other options keep their defaults, including the merging of surface
vertices closer than half the spacing.

The fine spacing makes the references nearly converged: from 0.25 to 0.05
bohr Multiwfn's water area changes by 1e-3 relative and its variance by
5e-4, and at 0.05 bohr the two programs agree to about 3e-5 on water's area,
variance, Pi and MPI.

**Test.** The test `multiwfn-surface-reference` runs
`electrostatic_potential_surface()` with its default spacing of 0.25 bohr
and compares the volume, areas, global extrema, averages, variances, nu,
Pi and MPI, and for water and acrolein the local extrema (count, value and
position). Multiwfn prints variances, Pi and MPI with more digits in
kcal/mol than in hartree, so those values are read in kcal/mol and divided
by Multiwfn's 627.51 kcal/mol per hartree. The differences are kumo's
discretization error at 0.25 bohr, which falls as the spacing squared, plus
Multiwfn's at 0.05 bohr. The largest observed differences, and the
tolerances, are:

| Quantity | Largest difference | Tolerance |
|----------|--------------------|-----------|
| Volume and areas | 1.6e-3 relative (O2 negative area) | 5e-3 relative |
| Averages, Pi, MPI | 1.6e-3 relative (O2 positive average) | 5e-3 relative + 2e-6 hartree/e |
| Variances | 3.4e-3 relative (O2) | 1e-2 relative + 2e-8 hartree^2 |
| nu | 2.2e-4 (O2) | 1e-3 |
| Nonpolar area fraction | 1.9e-4 (water) | 2e-3 |
| Extremum values | 1.4e-5 hartree/e (acrolein) | 4e-5 hartree/e |
| Extremum positions | 0.13 bohr (acrolein) | 0.4 bohr |

O2 has the largest relative differences because its potential is small
(below 6 kcal/mol in magnitude): kumo's own variance changes by 5e-3 from
0.25 to 0.05 bohr. The local extrema of O2 are not compared: its
potential is constant along rings around the axis, and Multiwfn lists 141
minima and 90 maxima on them, an artefact of the mesh.

## Multiwfn references for the average local ionization energy

The same analysis with the average local ionization energy (ALIE) as the
mapped function, on the same isosurface, with the same build and settings,
on 2026-11-13:

| File | FCHK |
|------|------|
| `h2o_alie_surface.txt` | `hybrids/h2o_b3lyp_def2tzvp.fchk.gz` |
| `o2_alie_surface.txt` | `o2.fchk.gz` |
| `acrolein_alie_surface.txt` | `acrolein.fchk.gz` |

The files are trimmed in the same way. The water analysis was run with this
input:

    h2o.fchk
    12
    2
    2
    3
    0.05
    0
    -1
    -1
    q

where `2` and `2` select the ALIE as the mapped function. Multiwfn's ALIE
(`avglocion` in `function.f90`) is sum_i n_i |e_i| phi_i^2 / sum_i n_i
phi_i^2 over the occupied orbitals, the definition kumo uses; for the ALIE
it prints the average and the variance about it, in hartree.

At 0.05 bohr kumo and Multiwfn agree to 4e-7 relative on the average of
water and acrolein, 2.5e-4 on their variances, 5e-4 on that of O2, and
5e-5 hartree on every local maximum of acrolein.

**Test.** The test `multiwfn-surface-reference` runs
`local_ionization_energy_surface()` with its default spacing of 0.25 bohr
and compares the volume, area, global extrema, average and variance, and
for water and acrolein the local minima (count, value and position):

| Quantity | Largest difference | Tolerance |
|----------|--------------------|-----------|
| Volume and area | as for the ESP | 5e-3 relative |
| Average | 1.5e-5 relative (acrolein) | 5e-5 relative |
| Variance | 6.1e-3 relative (O2) | 2e-2 relative |
| Global minimum | 1.0e-5 hartree (acrolein) | 4e-5 hartree |
| Global maximum | 1.4e-4 hartree (water) | 4e-4 hartree |
| Local minima, values | 2.5e-4 hartree (acrolein) | 8e-4 hartree |
| Local minima, positions | 0.14 bohr (acrolein) | 0.4 bohr |

The variance of O2 moves by 7e-3 in kumo from 0.25 to 0.05 bohr, and the
acrolein minimum with the largest difference converges to Multiwfn's value
within 1e-6 hartree at 0.05 bohr, so these are discretization errors.
Local maxima are not compared: they lie on sharp ridges at the hydrogen
tips, and at 0.25 bohr kumo finds 9 maxima for acrolein where Multiwfn
finds 8 at 0.05 bohr. As for the ESP, the local extrema of O2 lie on rings
and are not compared.

**Citation.** Multiwfn asks that work using it cite T. Lu and F. Chen,
J. Comput. Chem. 33, 580 (2012) and T. Lu, J. Chem. Phys. 161, 082503
(2024).
