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
if "$KUMO" "$FCHK" --spacing -1 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Invalid spacing unexpectedly accepted' >&2; exit 1
fi
if "$KUMO" "$FCHK" --mo 9999 --spacing 1 --output "$TMP/bad" > "$TMP/error.log" 2>&1; then
    echo 'Invalid orbital unexpectedly accepted' >&2; exit 1
fi
"$KUMO" --version | grep -q '1.3.0-dev'
echo 'CLI end-to-end OK'
