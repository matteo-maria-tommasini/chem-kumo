#!/bin/sh
# kumo --list-mo: the spin header names the expected wavefunction type, and
# for a restricted wavefunction the HOMO/LUMO labels sit on the right orbitals.
set -eu
fail() { echo "FAIL: $*" >&2; exit 1; }
[ $# -eq 3 ] || { echo "usage: $0 KUMO FILE.fchk restricted|alpha" >&2; exit 2; }
KUMO=$1
FCHK=$2
EXPECTED=$3
OUTPUT=$($KUMO "$FCHK" --list-mo) || fail "kumo $FCHK --list-mo exited with status $?"
printf '%s\n' "$OUTPUT" | grep -F "Spin: $EXPECTED" >/dev/null ||
    fail "--list-mo output does not contain 'Spin: $EXPECTED'"
if [ "$EXPECTED" = restricted ]; then
    printf '%s\n' "$OUTPUT" | grep -E '^ *15 +2 +.*HOMO$' >/dev/null ||
        fail "--list-mo does not mark orbital 15 (occupation 2) as HOMO"
    printf '%s\n' "$OUTPUT" | grep -E '^ *16 +0 +.*LUMO$' >/dev/null ||
        fail "--list-mo does not mark orbital 16 (occupation 0) as LUMO"
fi
