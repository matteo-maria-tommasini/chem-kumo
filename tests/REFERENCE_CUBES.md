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
