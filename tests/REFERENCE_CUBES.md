# Gaussian Cube reference files

The files `o2_density.cube` and `acrolein_density.cube` bundled with the
uploaded `chem-kumo-0.2.1-dev.tgz` are **not** the full-density references:
their numerical integrals are approximately 12.0046 and 22.0149 electrons.
The full-density references must be regenerated with:

    gcubegen09d 0 FDensity=SCF o2.fchk o2_density.cube -3 h
    gcubegen09d 0 FDensity=SCF acrolein.fchk acrolein_density.cube -3 h

and should integrate to approximately 16.0872 and 30.3105 electrons on
these grids. Replace the bundled cubes before expecting the tests to pass.
The test intentionally checks the reference integral against the FCHK
number of electrons to prevent accidental use of valence-only references.
