#!/bin/sh
set -eu
KUMO=$1
FCHK=$2
TMP=$(mktemp -d "${TMPDIR:-/tmp}/kumo-cli.XXXXXXXX")
trap 'rm -rf "$TMP"' EXIT HUP INT TERM
PREFIX="$TMP/o2"
"$KUMO" "$FCHK" --spacing 1.0 --margin 1.0 --output "$PREFIX" --spin all > "$TMP/density.log"
for kind in density alpha beta spin; do
    file="${PREFIX}_${kind}.cube"
    test -s "$file" || { echo "Missing density cube: $file" >&2; exit 1; }
    # Cube header: two comments, natoms/origin, three nonempty grid vectors.
    awk 'NR==3 {if (NF<4 || $1<=0) exit 1} NR>=4 && NR<=6 {if (NF<4 || $1<=0) exit 1} END {if (NR<=7) exit 1}' "$file" || exit 1
done
"$KUMO" "$FCHK" --mo frontier --spin both --spacing 1.0 --margin 1.0 --output "$PREFIX" > "$TMP/orbitals.log"
for name in alpha_mo_9 alpha_mo_10 beta_mo_7 beta_mo_8; do
    file="${PREFIX}_${name}.cube"
    test -s "$file" || { echo "Missing orbital cube: $file" >&2; exit 1; }
    awk 'NR==3 {if (NF<4 || $1>=0) exit 1} NR>=4 && NR<=6 {if (NF<4 || $1<=0) exit 1} END {if (NR<=8) exit 1}' "$file" || exit 1
done
# Density properties: --property density reproduces the default output
# exactly, and derivative Cube files have valid headers and plausible signs
# (|grad rho| and tau of the total density are nonnegative, the Laplacian of
# the total density is negative near the nuclei).
"$KUMO" "$FCHK" --spacing 1.0 --margin 1.0 --output "$TMP/default" > "$TMP/default.log"
"$KUMO" "$FCHK" --spacing 1.0 --margin 1.0 --output "$TMP/explicit" --property density > "$TMP/explicit.log"
tail -n +3 "$TMP/default_density.cube" > "$TMP/default.values"
tail -n +3 "$TMP/explicit_density.cube" > "$TMP/explicit.values"
cmp -s "$TMP/default.values" "$TMP/explicit.values" || {
    echo '--property density differs from the default density output' >&2; exit 1; }
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin all \
    --property gradient,laplacian,tau > "$TMP/derivatives.log"
test "$(wc -l < "$TMP/derivatives.log")" -eq 12 || {
    echo 'Expected 12 derivative Cube files' >&2; cat "$TMP/derivatives.log" >&2; exit 1; }
for kind in density alpha beta spin; do
    for property in gradient laplacian tau; do
        file="${PREFIX}_${kind}_${property}.cube"
        grep -Fqx "$file" "$TMP/derivatives.log" || { echo "Not reported: $file" >&2; exit 1; }
        test -s "$file" || { echo "Missing derivative cube: $file" >&2; exit 1; }
        awk 'NR==3 {if (NF<4 || $1<=0) exit 1} NR>=4 && NR<=6 {if (NF<4 || $1<=0) exit 1} END {if (NR<=7) exit 1}' "$file" || exit 1
    done
done
values() { awk 'NR == 3 {skip = 6 + ($1 < 0 ? -$1 : $1)} NR > 6 && NR > skip {for (i = 1; i <= NF; ++i) print $i}' "$1"; }
values "${PREFIX}_density_gradient.cube" | awk '$1 < 0 {exit 1}' || {
    echo 'Negative gradient norm' >&2; exit 1; }
values "${PREFIX}_density_tau.cube" | awk '$1 < 0 {exit 1}' || {
    echo 'Negative total tau' >&2; exit 1; }
values "${PREFIX}_density_laplacian.cube" | awk '$1 < 0 {found = 1} END {exit !found}' || {
    echo 'Total-density Laplacian has no negative values' >&2; exit 1; }
# Derived properties: ELF of the total, alpha and beta densities in [0, 1];
# the reduced gradient in [0, 100]; sign(lambda2) rho with negative values.
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin total \
    --property elf,rdg,sign-lambda2-rho > "$TMP/derived.log"
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin alpha \
    --property elf >> "$TMP/derived.log"
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --spin beta \
    --property elf >> "$TMP/derived.log"
for file in density_elf density_rdg density_sign-lambda2-rho alpha_elf beta_elf; do
    grep -Fqx "${PREFIX}_${file}.cube" "$TMP/derived.log" || {
        echo "Not reported: ${PREFIX}_${file}.cube" >&2; exit 1; }
done
for file in density_elf alpha_elf beta_elf; do
    values "${PREFIX}_${file}.cube" | awk '$1 < 0 || $1 > 1 {exit 1}' || {
        echo "$file outside [0, 1]" >&2; exit 1; }
done
values "${PREFIX}_density_rdg.cube" | awk '$1 < 0 || $1 > 100 {exit 1}' || {
    echo 'Reduced gradient outside [0, 100]' >&2; exit 1; }
values "${PREFIX}_density_sign-lambda2-rho.cube" | awk '$1 < 0 {found = 1} END {exit !found}' || {
    echo 'sign(lambda2) rho has no negative values' >&2; exit 1; }
# Electrostatic potential: this grid passes through both O2 nuclei, where the
# value is the potential at the nucleus (finite). The Cube comment names the
# method, whatever the density backend.
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$PREFIX" --property esp > "$TMP/esp.log"
grep -Fqx "${PREFIX}_density_esp.cube" "$TMP/esp.log" || {
    echo "Not reported: ${PREFIX}_density_esp.cube" >&2; exit 1; }
sed -n 2p "${PREFIX}_density_esp.cube" | grep -Fq 'Coulomb integrals over the AO density matrix, error <= 1e-08 hartree/e' || {
    echo 'ESP Cube comment does not name the method' >&2; exit 1; }
values "${PREFIX}_density_esp.cube" | awk '$1 + 0 != $1 || $1 > 1e3 {exit 1} $1 > 1 {found = 1} END {exit !found}' || {
    echo 'ESP values not finite, or none above 1 hartree/e near the nuclei' >&2; exit 1; }
# --esp-tolerance 0 is exact; the printed values agree with the default
# (error bound 1e-8) to the six digits of the Cube file, give or take one
# unit of the last digit.
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/exact" --property esp \
    --esp-tolerance 0 > "$TMP/esp-exact.log"
sed -n 2p "$TMP/exact_density_esp.cube" | grep -Fq 'exact Coulomb integrals' || {
    echo 'Exact ESP Cube comment does not say exact' >&2; exit 1; }
paste "$TMP/exact_density_esp.cube" "${PREFIX}_density_esp.cube" | awk 'NR > 8 {
    n = NF / 2; for (i = 1; i <= n; ++i) { d = $i - $(i + n); if (d < 0) d = -d;
    s = ($i < 0 ? -$i : $i); if (d > 1e-5 * s + 1e-12) exit 1 } }' || {
    echo 'ESP with --esp-tolerance 0 differs from the default' >&2; exit 1; }
if "$KUMO" "$FCHK" --spacing 1 --esp-tolerance 1e-6 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo '--esp-tolerance without --property esp unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --spacing 1 --property esp --esp-tolerance -1 --output "$TMP/bad" \
    > "$TMP/error.log" 2>&1; then
    echo 'Negative --esp-tolerance unexpectedly accepted' >&2; exit 1
fi
# Density backends: the Cube comment names the contraction used, and the AO-MO
# and density-matrix contractions agree (O2 is AO-MO eligible, so AUTO
# selects AO-MO).
for backend in auto ao-p ao-mo; do
    "$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/$backend" --backend "$backend" \
        --property density,laplacian > "$TMP/$backend.log"
done
for backend in auto ao-mo; do
    sed -n 2p "$TMP/${backend}_density.cube" | grep -Fq 'AO-MO contraction' || {
        echo "--backend $backend: Cube comment does not name the AO-MO contraction" >&2; exit 1; }
done
sed -n 2p "$TMP/ao-p_density_laplacian.cube" | grep -Fq 'AO density-matrix contraction' || {
    echo '--backend ao-p: Cube comment does not name the density-matrix contraction' >&2; exit 1; }
for name in density density_laplacian; do
    values "$TMP/ao-p_${name}.cube" > "$TMP/ao-p.values"
    values "$TMP/ao-mo_${name}.cube" > "$TMP/ao-mo.values"
    paste "$TMP/ao-p.values" "$TMP/ao-mo.values" | awk '
        { d = $1 - $2; if (d < 0) d = -d; s = ($1 < 0 ? -$1 : $1); if (s < 1) s = 1
          if (d > 1e-6 * s) bad = 1 }
        END { exit bad }' || { echo "AO-P and AO-MO $name Cube files disagree" >&2; exit 1; }
done
# --verbose: progress messages go to standard error only; standard output and
# the Cube files are the same as without it. Standard error is not a
# terminal here, so no progress bar (carriage returns) is drawn.
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/quiet" --spin all \
    --property density,tau > "$TMP/quiet.out" 2> "$TMP/quiet.err"
"$KUMO" "$FCHK" --spacing 0.5 --margin 1.0 --output "$TMP/loud" --spin all \
    --property density,tau --verbose > "$TMP/loud.out" 2> "$TMP/loud.err"
test ! -s "$TMP/quiet.err" || { echo 'Messages on stderr without --verbose' >&2; exit 1; }
sed "s|$TMP/quiet|$TMP/loud|" "$TMP/quiet.out" | cmp -s - "$TMP/loud.out" || {
    echo '--verbose changes standard output' >&2; exit 1; }
for message in 'Reading formatted checkpoint' 'basis functions' 'Grid:' \
    'Density contraction:' 'Evaluating the value and kinetic energy density tau' \
    'Evaluated in' 'Writing 8 Cube files' 'Finished in'; do
    grep -Fq "$message" "$TMP/loud.err" || {
        echo "--verbose message missing: $message" >&2; cat "$TMP/loud.err" >&2; exit 1; }
done
if tr -d '\n' < "$TMP/loud.err" | grep -q "$(printf '\r')"; then
    echo '--verbose drew a progress bar on a non-terminal' >&2; exit 1
fi
for name in density alpha beta spin density_tau alpha_tau beta_tau spin_tau; do
    cmp -s "$TMP/quiet_${name}.cube" "$TMP/loud_${name}.cube" || {
        echo "--verbose changes ${name} Cube file" >&2; exit 1; }
done
"$KUMO" "$FCHK" --mo homo -v --spacing 1.0 --margin 1.0 --output "$TMP/mo" \
    > "$TMP/mo.out" 2> "$TMP/mo.err"
grep -Fq 'Evaluating alpha molecular orbital 9' "$TMP/mo.err" || {
    echo '--verbose --mo message missing' >&2; cat "$TMP/mo.err" >&2; exit 1; }

if "$KUMO" "$FCHK" --spacing 1 --backend mo --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Unknown --backend unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --mo homo --backend ao-p --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo '--backend with --mo unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --spacing 1 --property hessian --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Unknown --property unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --mo homo --property tau --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo '--property with --mo unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --spacing 1 --spin difference --property elf --output "$TMP/bad" \
    > "$TMP/error.log" 2>&1; then
    echo 'ELF of the spin difference unexpectedly accepted' >&2; exit 1
fi
for property in rdg sign-lambda2-rho esp; do
    if "$KUMO" "$FCHK" --spacing 1 --spin alpha --property "$property" --output "$TMP/bad" \
        > "$TMP/error.log" 2>&1; then
        echo "--property $property with --spin alpha unexpectedly accepted" >&2; exit 1
    fi
done
if "$KUMO" "$FCHK" --spacing 1 --property tau, --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Trailing comma in --property unexpectedly accepted' >&2; exit 1
fi

if "$KUMO" "$FCHK" --spacing -1 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Invalid spacing unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --mo 9999 --spacing 1 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Invalid orbital unexpectedly accepted' >&2; exit 1
fi

EXPECTED_VERSION=$3
ACTUAL_VERSION=$("$KUMO" --version)

printf '%s\n' "$ACTUAL_VERSION" |
    grep -Fq -- "$EXPECTED_VERSION" || {
        echo "Unexpected CLI version: $ACTUAL_VERSION" >&2
        echo "Expected version: $EXPECTED_VERSION" >&2
        exit 1
    }

echo 'CLI end-to-end OK'
