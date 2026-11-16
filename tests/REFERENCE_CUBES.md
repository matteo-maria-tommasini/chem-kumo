# Gaussian Cube reference files

This directory contains Gaussian 09 reference Cube files used for
quantitative validation of chem-kumo.

The electron-density reference files were generated using the
full SCF electron density:

    cubegen 0 FDensity=SCF o2.fchk o2_density.cube -3 h
    cubegen 0 FDensity=SCF acrolein.fchk acrolein_density.cube -3 h

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

**Citation.** The Gaussian FCHK files in `tests/data/` and the Gaussian
reference Cube files were computed with Gaussian 09, Revision D.01 (the Cube
files with its `cubegen` utility). Work that uses these data should cite:

> Gaussian 09, Revision D.01,
> M. J. Frisch, G. W. Trucks, H. B. Schlegel, G. E. Scuseria,
> M. A. Robb, J. R. Cheeseman, G. Scalmani, V. Barone, B. Mennucci,
> G. A. Petersson, H. Nakatsuji, M. Caricato, X. Li, H. P. Hratchian,
> A. F. Izmaylov, J. Bloino, G. Zheng, J. L. Sonnenberg, M. Hada,
> M. Ehara, K. Toyota, R. Fukuda, J. Hasegawa, M. Ishida, T. Nakajima,
> Y. Honda, O. Kitao, H. Nakai, T. Vreven, J. A. Montgomery, Jr.,
> J. E. Peralta, F. Ogliaro, M. Bearpark, J. J. Heyd, E. Brothers,
> K. N. Kudin, V. N. Staroverov, T. Keith, R. Kobayashi, J. Normand,
> K. Raghavachari, A. Rendell, J. C. Burant, S. S. Iyengar, J. Tomasi,
> M. Cossi, N. Rega, J. M. Millam, M. Klene, J. E. Knox, J. B. Cross,
> V. Bakken, C. Adamo, J. Jaramillo, R. Gomperts, R. E. Stratmann,
> O. Yazyev, A. J. Austin, R. Cammi, C. Pomelli, J. W. Ochterski,
> R. L. Martin, K. Morokuma, V. G. Zakrzewski, G. A. Voth,
> P. Salvador, J. J. Dannenberg, S. Dapprich, A. D. Daniels,
> O. Farkas, J. B. Foresman, J. V. Ortiz, J. Cioslowski,
> and D. J. Fox, Gaussian, Inc., Wallingford CT, 2013.

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

### Orbital-energy fields of 1.13.0

Cube files of the ALIE and of Multiwfn's user functions 27 (local electron
affinity), -27 (local electron attachment energy), 28 (local
electronegativity) and 29 (local hardness) were computed on 2026-11-13 with
the same Multiwfn build, grid extension (3 bohr), spacing (0.3 bohr), input
sequence and scratch `settings.ini` (with `ipolarpara= 0` and `iuserfunc`
set to the function), from `acrolein.fchk.gz` (restricted; one virtual
orbital with negative energy, the LUMO at -0.0115 hartree) and `o2.fchk.gz`
(unrestricted; two beta virtual orbitals at -0.0464 hartree).

| File | FCHK | Multiwfn function |
|------|------|-------------------|
| `acrolein_alie.cube.gz` | acrolein | 18, ALIE |
| `acrolein_lea.cube.gz` | acrolein | 100 with `iuserfunc= 27` |
| `acrolein_leae.cube.gz` | acrolein | 100 with `iuserfunc= -27` |
| `acrolein_local_electronegativity.cube.gz` | acrolein | 100 with `iuserfunc= 28` |
| `acrolein_local_hardness.cube.gz` | acrolein | 100 with `iuserfunc= 29` |
| `o2_alie.cube.gz` | O2 | 18, ALIE |
| `o2_lea.cube.gz` | O2 | 100 with `iuserfunc= 27` |
| `o2_leae.cube.gz` | O2 | 100 with `iuserfunc= -27` |

The test `multiwfn-reference` compares them point by point:

| Quantity | Tolerance | Largest difference observed |
|----------|-----------|-----------------------------|
| ALIE, LEA, local hardness | 2e-5 relative | 5.0e-6 relative |
| LEAE | 2e-5 relative + 1e-10 hartree; 2e-5 relative alone where \|LEAE\| >= 1e-3 hartree | 6.4e-6 relative where \|LEAE\| >= 1e-3 (12258 acrolein and 11788 O2 points); 5.4e-11 hartree near the nodal planes of the virtual orbitals |
| Local electronegativity | 2e-5 of (\|ALIE\| + \|LEA\|) / 2 | 2.7e-6 of that scale, where the electronegativity crosses zero |

LEAE vanishes on the nodal planes of the virtual orbitals with negative
energy. There both programs give rounding noise (-4e-39 versus -1e-34), and
close to them, where LEAE is a few 1e-6 hartree, the relative precision of
the cancelling orbital sums drops in both: the absolute term of 1e-10
hartree (1e-8 of the largest |LEAE|) covers these points. The test also
requires at least 10000 points with |LEAE| >= 1e-3 hartree and the
relative tolerance alone there, so it cannot pass on near-zero values only.

The largest difference is 0.54 of its tolerance (O2 LEAE).

### Energy densities of 2.1.0

Cube files of the Lagrangian kinetic energy density G (Multiwfn function 7,
kumo's tau) and of Multiwfn's user functions 81-83 (K_x, K_y, K_z), 84-86
(G_x, G_y, G_z), 35 (|V| / G), 17 (H / rho) and 13 (G / rho) were computed
on 2026-11-16 with the same Multiwfn build (`Multiwfn_noGUI`, 2026.9.1
macOS port), grid extension (3 bohr), spacing (0.3 bohr) and input
sequence, from a scratch copy of the distributed `settings.ini` with only
`iuserfunc` set to the function.

| File | FCHK | Multiwfn function |
|------|------|-------------------|
| `h2o_kinetic_g.cube.gz` | water | 7, G(r) |
| `o2_kinetic_g.cube.gz` | O2 | 7, G(r) |
| `h2o_kinetic_g_x.cube.gz`, `_y`, `_z` | water | 100 with `iuserfunc= 84`, `85`, `86` |
| `h2o_kinetic_k_x.cube.gz`, `_y`, `_z` | water | 100 with `iuserfunc= 81`, `82`, `83` |
| `o2_kinetic_k_x.cube.gz` | O2 | 100 with `iuserfunc= 81` |
| `h2o_potential_kinetic_ratio.cube.gz` | water | 100 with `iuserfunc= 35`, \|V\| / G |
| `o2_potential_kinetic_ratio.cube.gz` | O2 | 100 with `iuserfunc= 35` |
| `h2o_energy_density_per_electron.cube.gz` | water | 100 with `iuserfunc= 17`, H / rho |
| `h2o_kinetic_g_per_electron.cube.gz` | water | 100 with `iuserfunc= 13`, G / rho |

Multiwfn computes G_a as (1/2) sum_i n_i (d_a psi_i)^2, K_a as
-(1/2) sum_i n_i psi_i d_a^2 psi_i, |V| / G as |-K - G| / G, H / rho as
-K / rho and G / rho directly (`Lagkin` and `Hamkin` in `function.f90`).
The test `multiwfn-reference` compares them like the other files:

| Field | Tolerance | Largest difference |
|-------|-----------|--------------------|
| G (water, O2), G / rho | 2e-5 relative | 7.4e-6 relative |
| G_x, G_y, G_z | 5e-5 relative | 1.7e-5 (water G_x, where it is 0.02) |
| K_x, K_y, K_z (water, O2 K_x) | 5e-5 of \|G_a\| + \|H_aa\| / 4 | 1.2e-5 of that scale |
| \|V\| / G (water, O2) | 5e-5 of (2 \|G\| + \|lap rho\| / 4) / G | 2.4e-5 of that scale (O2) |
| H / rho | 5e-5 of (\|G\| + \|lap rho\| / 4) / rho | 9.0e-6 of that scale |

H_aa is the diagonal of the density Hessian. K_a, |V| / G and H / rho
cancel terms as K, V and H do, and are compared on the scale of those
terms with the tolerance of K, V and H. The O2 grid contains one nucleus,
at (0, 0, -1.143284): the six decimals of the Cube header move it by up to
5e-7 bohr, and on the nuclear cusp G changes by about 2Z = 32 relative per
bohr, so G differs there by 2.9e-5. That point is skipped for O2 G (points
within 1e-3 bohr of a nucleus, as for the ESP); every other point agrees
within 5.5e-6. The largest difference is 0.47 of its tolerance
(O2 |V| / G).

### Orbital-weighted Fukui functions of 2.1.0

Cube files of Multiwfn's user functions 95 (orbital-weighted f+), 96 (f-),
97 (f0) and 98 (dual descriptor) were computed on 2026-11-16 with the same
Multiwfn build (`Multiwfn_noGUI`, 2026.9.1 macOS port), grid extension (3 bohr), spacing (0.3 bohr), input sequence
(function 100) and scratch `settings.ini` (with `ipolarpara= 0` and
`iuserfunc` set to the function), with Multiwfn's default width
delta = 0.1 hartree, from `acrolein.fchk.gz` and
`hybrids/h2o_b3lyp_def2tzvp.fchk.gz` (both restricted closed-shell;
Multiwfn defines these functions only for closed shells).

| File | FCHK | Multiwfn function |
|------|------|-------------------|
| `acrolein_ow_fukui_plus.cube.gz` | acrolein | 100 with `iuserfunc= 95` |
| `acrolein_ow_fukui_minus.cube.gz` | acrolein | 100 with `iuserfunc= 96` |
| `acrolein_ow_fukui_zero.cube.gz` | acrolein | 100 with `iuserfunc= 97` |
| `acrolein_ow_dual_descriptor.cube.gz` | acrolein | 100 with `iuserfunc= 98` |
| `h2o_ow_fukui_plus.cube.gz` | water | 100 with `iuserfunc= 95` |
| `h2o_ow_fukui_minus.cube.gz` | water | 100 with `iuserfunc= 96` |

Multiwfn (`orbwei_Fukui` in `function.f90`) keeps only the orbitals whose
unnormalized weight exp(-((mu - e_i) / delta)^2) is at least 1e-3,
scanning up from the LUMO (f+) and down from the HOMO (f-) and stopping at
the first orbital below; kumo sums over all orbitals (docs/API.md). At
delta = 0.1 Multiwfn keeps the LUMO alone and the HOMO and HOMO-1 for
acrolein, which drops 0.44% of the weight of f+ and 0.30% of f-, and the
LUMO, LUMO+1, HOMO and HOMO-1 for water (9e-11 and 5.7e-6 of the weight).
These differences are far above the Cube precision, so the test
`multiwfn-reference` rebuilds Multiwfn's truncated sums from the orbitals
evaluated by kumo, adds the dropped part (kumo's f+ minus the truncated
sum, and likewise for f-) to each Cube value, and compares the result with
kumo's field:

| Quantity | Tolerance | Largest difference observed |
|----------|-----------|-----------------------------|
| f+, f-, f0 | 2e-5 relative + 1e-8 bohr^-3 | 9.4e-6 relative where \|value\| >= 1e-4 (water f+); 3e-9 bohr^-3 below |
| Dual descriptor | 2e-5 of \|f+\| + \|f-\|, + 1e-8 bohr^-3 | 5.7e-6 of that scale where \|value\| >= 1e-4 |

Where |value| >= 1e-4 bohr^-3 (6659 to 15591 points per Cube; the largest
values are 0.14 to 0.47) the relative tolerance must hold alone. Below, at
seven points with values of 1e-6 to 7e-5, the two programs differ by up to
3e-9 bohr^-3 (4e-5 relative); the absolute term, about three times that,
covers them. The correction matters: at some points (on the nodal surfaces
of the kept orbitals) it is the whole of kumo's value, which Multiwfn
gives as nearly zero. The largest difference is 0.47 of its tolerance
(0.28 overall, 0.47 where the relative tolerance holds alone).

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

### LEAE and LEA mapped on the surface (1.13.0)

The same analysis with the local electron attachment energy (mapped
function -4) and the local electron affinity (mapped function 4), on
2026-11-13, with the input sequence of the ALIE (`12`, `2`, then `-4` or
`4`, `3`, `0.05`, `0`, `-1`, `-1`, `q`):

| File | FCHK | Function | Isosurface |
|------|------|----------|------------|
| `acrolein_leae_surface.txt` | `acrolein.fchk.gz` | -4, LEAE | 0.004 |
| `o2_leae_surface.txt` | `o2.fchk.gz` | -4, LEAE | 0.004 |
| `acrolein_lea_surface.txt` | `acrolein.fchk.gz` | 4, LEA | 0.001 |

Choosing LEAE makes Multiwfn change the isovalue to 0.004 ("Isovalue has
been changed to 0.004 a.u."), the surface of the original work, so the test
maps LEAE on rho = 0.004. The files are trimmed as above. The LEAE maxima
(368 for acrolein) are the zeros on the nodal plane of the pi* orbital and
are not compared; the extrema of O2 lie on rings.

The test `multiwfn-surface-reference` maps them with
`map_on_density_isosurface()` at 0.25 bohr:

| Quantity | Largest difference | Tolerance |
|----------|--------------------|-----------|
| Volume and area | 6.3e-4 relative (O2) | 5e-3 relative |
| LEAE average | 3.0e-4 relative (acrolein) | 1e-3 relative |
| LEAE variance | 7.6e-3 relative (O2) | 2.5e-2 relative |
| LEAE global minimum | 4.6e-7 hartree | 1.5e-6 hartree |
| LEAE local minima (acrolein) | 1.8e-5 hartree, 0.07 bohr | 6e-5 hartree, 0.4 bohr |
| LEA average | 5.1e-5 relative | 1.5e-4 relative |
| LEA variance | 1.7e-3 relative | 5e-3 relative |
| LEA global extrema | 1.2e-5 hartree | 4e-5 hartree |
| LEA local minima | 1.9e-3 hartree, 0.15 bohr | 6e-3 hartree, 0.4 bohr |

At 0.05 bohr kumo's LEAE variances agree with Multiwfn's to 1e-3 and its
LEA minima to 3e-5 hartree, so the differences at 0.25 bohr are
discretization errors; the LEA minima lie on ridges near the hydrogens, as
the ALIE maxima. The largest difference is 0.36 of its tolerance.

**Citation.** Multiwfn asks that work using it cite T. Lu and F. Chen,
J. Comput. Chem. 33, 580 (2012) and T. Lu, J. Chem. Phys. 161, 082503
(2024).
