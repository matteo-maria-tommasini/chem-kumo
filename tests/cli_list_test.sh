#!/bin/sh
set -eu
KUMO=$1
FCHK=$2
EXPECTED=$3
OUTPUT=$($KUMO "$FCHK" --list-mo)
printf '%s\n' "$OUTPUT" | grep -F "Spin: $EXPECTED" >/dev/null
if [ "$EXPECTED" = restricted ]; then
    printf '%s\n' "$OUTPUT" | grep -E '^ *15 +2 +.*HOMO$' >/dev/null
    printf '%s\n' "$OUTPUT" | grep -E '^ *16 +0 +.*LUMO$' >/dev/null
fi
