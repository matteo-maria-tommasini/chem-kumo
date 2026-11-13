#!/bin/sh
set -eu
[ $# -eq 3 ] || { echo "usage: $0 KUMO FILE.fchk EXPECTED_VERSION" >&2; exit 2; }
KUMO=$1
FCHK=$2
TMP=$(mktemp -d "${TMPDIR:-/tmp}/kumo-cli.XXXXXXXX")
FAILED=
# Every failure prints a message: through fail(), or, for a command that
# fails under set -e, through the EXIT trap.
fail() { FAILED=1; echo "FAIL: $*" >&2; exit 1; }
trap 'status=$?; rm -rf "$TMP"
      if [ $status -ne 0 ] && [ -z "$FAILED" ]; then
          echo "FAIL: unexpected error (exit status $status)" >&2; fi' EXIT
trap 'exit 1' HUP INT TERM
# Run kumo; a nonzero exit status is a failure.
run_kumo() { "$KUMO" "$@" || fail "kumo $* exited with status $?"; }
# Cube header: two comments, the atom count and origin (atom count > 0 for a
# density, < 0 for an orbital Cube, which has one more header line), three
# nonempty grid vectors, then data.
check_cube_header() { # FILE density|orbital
    awk -v orbital="$([ "$2" = orbital ] && echo 1 || echo 0)" '
        NR==3 {if (NF<4 || (orbital ? $1>=0 : $1<=0)) exit 1}
        NR>=4 && NR<=6 {if (NF<4 || $1<=0) exit 1}
        END {if (NR<=(orbital ? 8 : 7)) exit 1}' "$1" || fail "invalid $2 Cube header: $1"
}
PREFIX="$TMP/o2"
run_kumo "$FCHK" --spacing 1.0 --margin 1.0 --output "$PREFIX" --spin all > "$TMP/density.log"
for kind in density alpha beta spin; do
    file="${PREFIX}_${kind}.cube"
    test -s "$file" || fail "Missing density cube: $file"
    check_cube_header "$file" density
done
run_kumo "$FCHK" --mo frontier --spin both --spacing 1.0 --margin 1.0 --output "$PREFIX" > "$TMP/orbitals.log"
for name in alpha_mo_9 alpha_mo_10 beta_mo_7 beta_mo_8; do
    file="${PREFIX}_${name}.cube"
    test -s "$file" || fail "Missing orbital cube: $file"
    check_cube_header "$file" orbital
done
# Density properties: --property density reproduces the default output
# exactly, and derivative Cube files have valid headers and plausible signs
# (|grad rho| and tau of the total density are nonnegative, the Laplacian of
# the total density is negative near the nuclei).
run_kumo "$FCHK" --spacing 1.0 --margin 1.0 --output "$TMP/default" > "$TMP/default.log"
run_kumo "$FCHK" --spacing 1.0 --margin 1.0 --output "$TMP/explicit" --property density > "$TMP/explicit.log"
tail -n +3 "$TMP/default_density.cube" > "$TMP/default.values"
tail -n +3 "$TMP/explicit_density.cube" > "$TMP/explicit.values"
cmp -s "$TMP/default.values" "$TMP/explicit.values" ||
    fail '--property density differs from the default density output'
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin all \
    --property gradient,laplacian,tau > "$TMP/derivatives.log"
test "$(wc -l < "$TMP/derivatives.log")" -eq 12 || {
    cat "$TMP/derivatives.log" >&2; fail 'Expected 12 derivative Cube files'; }
for kind in density alpha beta spin; do
    for property in gradient laplacian tau; do
        file="${PREFIX}_${kind}_${property}.cube"
        grep -Fqx "$file" "$TMP/derivatives.log" || fail "Not reported: $file"
        test -s "$file" || fail "Missing derivative cube: $file"
        check_cube_header "$file" density
    done
done
values() { awk 'NR == 3 {skip = 6 + ($1 < 0 ? -$1 : $1)} NR > 6 && NR > skip {for (i = 1; i <= NF; ++i) print $i}' "$1"; }
values "${PREFIX}_density_gradient.cube" | awk '$1 < 0 {exit 1}' ||
    fail 'Negative gradient norm'
values "${PREFIX}_density_tau.cube" | awk '$1 < 0 {exit 1}' ||
    fail 'Negative total tau'
values "${PREFIX}_density_laplacian.cube" | awk '$1 < 0 {found = 1} END {exit !found}' ||
    fail 'Total-density Laplacian has no negative values'
# Derived properties: ELF of the total, alpha and beta densities in [0, 1];
# the reduced gradient in [0, 100]; sign(lambda2) rho with negative values.
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin total \
    --property elf,rdg,sign-lambda2-rho > "$TMP/derived.log"
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin alpha \
    --property elf >> "$TMP/derived.log"
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin beta \
    --property elf >> "$TMP/derived.log"
for file in density_elf density_rdg density_sign-lambda2-rho alpha_elf beta_elf; do
    grep -Fqx "${PREFIX}_${file}.cube" "$TMP/derived.log" ||
        fail "Not reported: ${PREFIX}_${file}.cube"
done
for file in density_elf alpha_elf beta_elf; do
    values "${PREFIX}_${file}.cube" | awk '$1 < 0 || $1 > 1 {exit 1}' ||
        fail "$file outside [0, 1]"
done
values "${PREFIX}_density_rdg.cube" | awk '$1 < 0 || $1 > 100 {exit 1}' ||
    fail 'Reduced gradient outside [0, 100]'
values "${PREFIX}_density_sign-lambda2-rho.cube" | awk '$1 < 0 {found = 1} END {exit !found}' ||
    fail 'sign(lambda2) rho has no negative values'
# Electrostatic potential: this grid passes through both O2 nuclei, where the
# value is the potential at the nucleus (finite). The Cube comment names the
# method, whatever the density backend.
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --property esp > "$TMP/esp.log"
grep -Fqx "${PREFIX}_density_esp.cube" "$TMP/esp.log" ||
    fail "Not reported: ${PREFIX}_density_esp.cube"
sed -n 2p "${PREFIX}_density_esp.cube" | grep -Fq 'Coulomb integrals over the AO density matrix, error <= 1e-08 hartree/e' ||
    fail 'ESP Cube comment does not name the method'
values "${PREFIX}_density_esp.cube" | awk '$1 + 0 != $1 || $1 > 1e3 {exit 1} $1 > 1 {found = 1} END {exit !found}' ||
    fail 'ESP values not finite, or none above 1 hartree/e near the nuclei'
# --esp-tolerance 0 is exact; the printed values agree with the default
# (error bound 1e-8) to the six digits of the Cube file, give or take one
# unit of the last digit.
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/exact" --property esp \
    --esp-tolerance 0 > "$TMP/esp-exact.log"
sed -n 2p "$TMP/exact_density_esp.cube" | grep -Fq 'exact Coulomb integrals' ||
    fail 'Exact ESP Cube comment does not say exact'
paste "$TMP/exact_density_esp.cube" "${PREFIX}_density_esp.cube" | awk 'NR > 8 {
    n = NF / 2; for (i = 1; i <= n; ++i) { d = $i - $(i + n); if (d < 0) d = -d;
    s = ($i < 0 ? -$i : $i); if (d > 1e-5 * s + 1e-12) exit 1 } }' ||
    fail 'ESP with --esp-tolerance 0 differs from the default'
if "$KUMO" "$FCHK" --spacing 1 --esp-tolerance 1e-6 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail '--esp-tolerance without --property esp unexpectedly accepted'
fi
if "$KUMO" "$FCHK" --spacing 1 --property esp --esp-tolerance -1 --output "$TMP/bad" \
    > "$TMP/error.log" 2>&1; then
    fail 'Negative --esp-tolerance unexpectedly accepted'
fi
# --esp-surface prints the potential statistics on the rho = 0.001 isosurface
# (O2: maximum about 5.63 kcal/mol); --isovalue, --spacing and --esp-tolerance
# change the surface and the potential.
run_kumo "$FCHK" --esp-surface > "$TMP/surface.log"
grep -Fqx 'Electrostatic potential on the total-density isosurface rho = 0.001 e/bohr^3' \
    "$TMP/surface.log" || fail '--esp-surface: no isosurface line'
grep -Eqx 'Surface: [0-9]+ vertices, [0-9]+ triangles, grid spacing 0.25 bohr' "$TMP/surface.log" ||
    fail '--esp-surface: no surface summary with the default spacing 0.25'
grep -Fq 'error <= 1e-08 hartree/e' "$TMP/surface.log" ||
    fail '--esp-surface: no potential error bound'
awk '$1 == "Maximum" && $3 > 5.6 && $3 < 5.7 {found = 1} END {exit !found}' "$TMP/surface.log" ||
    fail '--esp-surface: maximum of O2 not near 5.63 kcal/mol'
grep -Eq '^Surface maxima: [0-9]+$' "$TMP/surface.log" || fail '--esp-surface: no maxima list'
run_kumo "$FCHK" --esp-surface --isovalue 0.002 --spacing 0.3 --esp-tolerance 0 \
    > "$TMP/surface2.log"
grep -Fq 'rho = 0.002 e/bohr^3' "$TMP/surface2.log" && grep -Fq 'grid spacing 0.3 bohr' "$TMP/surface2.log" &&
    grep -Fq 'density matrix, exact' "$TMP/surface2.log" ||
    fail '--esp-surface: --isovalue, --spacing or --esp-tolerance 0 not applied'
area() { awk '$1 == "Area" {print $2}' "$1"; }
awk -v a="$(area "$TMP/surface.log")" -v b="$(area "$TMP/surface2.log")" 'BEGIN {exit !(b < a)}' ||
    fail '--esp-surface: the rho = 0.002 surface is not smaller than the rho = 0.001 one'
for options in '--property esp' '--output x' '--margin 2' '--spin total' '--mo homo' \
    '--isovalue 1e-9' '--isovalue 0.001 --isovalue 0.002'; do
    # shellcheck disable=SC2086 # options are split on purpose
    if "$KUMO" "$FCHK" --esp-surface $options > "$TMP/error.log" 2>&1; then
        fail "--esp-surface $options unexpectedly accepted"
    fi
done
if "$KUMO" "$FCHK" --spacing 1 --isovalue 0.002 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail '--isovalue without --esp-surface unexpectedly accepted'
fi
# --alie-surface prints the average local ionization energy statistics on the
# same isosurface (O2: minimum about 15.40 eV, average about 16.32 eV).
run_kumo "$FCHK" --alie-surface > "$TMP/alie.log"
grep -Fqx 'Average local ionization energy on the total-density isosurface rho = 0.001 e/bohr^3' \
    "$TMP/alie.log" || fail '--alie-surface: no isosurface line'
grep -Eqx 'Surface: [0-9]+ vertices, [0-9]+ triangles, grid spacing 0.25 bohr' "$TMP/alie.log" ||
    fail '--alie-surface: no surface summary with the default spacing 0.25'
awk '$1 == "Minimum" && $3 > 15.3 && $3 < 15.5 {found = 1} END {exit !found}' "$TMP/alie.log" ||
    fail '--alie-surface: minimum of O2 not near 15.40 eV'
awk '$1 == "Average" && $3 > 16.2 && $3 < 16.4 {found = 1} END {exit !found}' "$TMP/alie.log" ||
    fail '--alie-surface: average of O2 not near 16.32 eV'
grep -Eq '^Surface minima: [0-9]+$' "$TMP/alie.log" || fail '--alie-surface: no minima list'
run_kumo "$FCHK" --alie-surface --isovalue 0.002 --spacing 0.3 > "$TMP/alie2.log"
grep -Fq 'rho = 0.002 e/bohr^3' "$TMP/alie2.log" && grep -Fq 'grid spacing 0.3 bohr' "$TMP/alie2.log" ||
    fail '--alie-surface: --isovalue or --spacing not applied'
awk -v a="$(area "$TMP/alie.log")" -v b="$(area "$TMP/alie2.log")" 'BEGIN {exit !(b < a)}' ||
    fail '--alie-surface: the rho = 0.002 surface is not smaller than the rho = 0.001 one'
for options in '--esp-surface' '--esp-tolerance 0' '--property esp' '--output x' '--margin 2' \
    '--spin total' '--mo homo' '--backend ao-p' '--isovalue 1e-9'; do
    # shellcheck disable=SC2086 # options are split on purpose
    if "$KUMO" "$FCHK" --alie-surface $options > "$TMP/error.log" 2>&1; then
        fail "--alie-surface $options unexpectedly accepted"
    fi
done
# Density backends: the Cube comment names the contraction used, and the AO-MO
# and density-matrix contractions agree (O2 is AO-MO eligible, so AUTO
# selects AO-MO).
for backend in auto ao-p ao-mo; do
    run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/$backend" --backend "$backend" \
        --property density,laplacian > "$TMP/$backend.log"
done
for backend in auto ao-mo; do
    sed -n 2p "$TMP/${backend}_density.cube" | grep -Fq 'AO-MO contraction' ||
        fail "--backend $backend: Cube comment does not name the AO-MO contraction"
done
sed -n 2p "$TMP/ao-p_density_laplacian.cube" | grep -Fq 'AO density-matrix contraction' ||
    fail '--backend ao-p: Cube comment does not name the density-matrix contraction'
for name in density density_laplacian; do
    values "$TMP/ao-p_${name}.cube" > "$TMP/ao-p.values"
    values "$TMP/ao-mo_${name}.cube" > "$TMP/ao-mo.values"
    paste "$TMP/ao-p.values" "$TMP/ao-mo.values" | awk '
        { d = $1 - $2; if (d < 0) d = -d; s = ($1 < 0 ? -$1 : $1); if (s < 1) s = 1
          if (d > 1e-6 * s) bad = 1 }
        END { exit bad }' || fail "AO-P and AO-MO $name Cube files disagree"
done
# --verbose: progress messages go to standard error only; standard output and
# the Cube files are the same as without it. Standard error is not a
# terminal here, so no progress bar (carriage returns) is drawn.
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/quiet" --spin all \
    --property density,tau > "$TMP/quiet.out" 2> "$TMP/quiet.err"
run_kumo "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/loud" --spin all \
    --property density,tau --verbose > "$TMP/loud.out" 2> "$TMP/loud.err"
test ! -s "$TMP/quiet.err" || fail 'Messages on stderr without --verbose'
sed "s|$TMP/quiet|$TMP/loud|" "$TMP/quiet.out" | cmp -s - "$TMP/loud.out" ||
    fail '--verbose changes standard output'
for message in 'Reading formatted checkpoint' 'basis functions' 'Grid:' \
    'Occupied orbitals vs FCHK density matrices: error relative to the total matrix' \
    'AO-MO needs both <= 1e-08' \
    'Density contraction:' 'Evaluating the value and kinetic energy density tau' \
    'Evaluated in' 'Writing 8 Cube files' 'Finished in'; do
    grep -Fq "$message" "$TMP/loud.err" || {
        cat "$TMP/loud.err" >&2; fail "--verbose message missing: $message"; }
done
if tr -d '\n' < "$TMP/loud.err" | grep -q "$(printf '\r')"; then
    fail '--verbose drew a progress bar on a non-terminal'
fi
for name in density alpha beta spin density_tau alpha_tau beta_tau spin_tau; do
    cmp -s "$TMP/quiet_${name}.cube" "$TMP/loud_${name}.cube" ||
        fail "--verbose changes ${name} Cube file"
done
grep -Fq '(not met)' "$TMP/loud.err" && fail '--verbose reports a failed orbital check for O2'
# A density matrix that the occupied orbitals do not reproduce: --verbose
# reports the failed check and the AO-P contraction.
gzip -dc "$FCHK" | awk '
    done == 0 && previous ~ /^Total SCF Density/ { $2 = sprintf("%.8E", $2 + 0.05); done = 1 }
    { print; previous = $0 }' > "$TMP/inconsistent.fchk"
run_kumo "$TMP/inconsistent.fchk" --spacing 1.0 --margin 1.0 -v --output "$TMP/inconsistent" \
    > /dev/null 2> "$TMP/inconsistent.err"
grep -Fq 'AO-MO needs both <= 1e-08 (not met)' "$TMP/inconsistent.err" &&
    grep -Fq 'the occupied orbitals do not reproduce' "$TMP/inconsistent.err" || {
    cat "$TMP/inconsistent.err" >&2; fail '--verbose does not report the failed orbital check'; }
run_kumo "$FCHK" --mo homo -v --spacing 1.0 --margin 1.0 --output "$TMP/mo" \
    > "$TMP/mo.out" 2> "$TMP/mo.err"
grep -Fq 'Evaluating alpha molecular orbital 9' "$TMP/mo.err" || {
    cat "$TMP/mo.err" >&2; fail '--verbose --mo message missing'; }

if "$KUMO" "$FCHK" --spacing 1 --backend mo --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail 'Unknown --backend unexpectedly accepted'
fi
if "$KUMO" "$FCHK" --mo homo --backend ao-p --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail '--backend with --mo unexpectedly accepted'
fi
if "$KUMO" "$FCHK" --spacing 1 --property hessian --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail 'Unknown --property unexpectedly accepted'
fi
if "$KUMO" "$FCHK" --mo homo --property tau --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail '--property with --mo unexpectedly accepted'
fi
if "$KUMO" "$FCHK" --spacing 1 --spin difference --property elf --output "$TMP/bad" \
    > "$TMP/error.log" 2>&1; then
    fail 'ELF of the spin difference unexpectedly accepted'
fi
for property in rdg sign-lambda2-rho esp; do
    if "$KUMO" "$FCHK" --spacing 1 --spin alpha --property "$property" --output "$TMP/bad" \
        > "$TMP/error.log" 2>&1; then
        fail "--property $property with --spin alpha unexpectedly accepted"
    fi
done
if "$KUMO" "$FCHK" --spacing 1 --property tau, --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail 'Trailing comma in --property unexpectedly accepted'
fi

if "$KUMO" "$FCHK" --spacing -1 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail 'Invalid spacing unexpectedly accepted'
fi
if "$KUMO" "$FCHK" --mo 9999 --spacing 1 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    fail 'Invalid orbital unexpectedly accepted'
fi

EXPECTED_VERSION=$3
ACTUAL_VERSION=$("$KUMO" --version)

printf '%s\n' "$ACTUAL_VERSION" |
    grep -Fq -- "$EXPECTED_VERSION" ||
    fail "unexpected CLI version: $ACTUAL_VERSION (expected $EXPECTED_VERSION)"

echo 'CLI end-to-end OK'
