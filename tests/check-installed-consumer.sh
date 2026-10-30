#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
TMP=$(mktemp -d "${TMPDIR:-/tmp}/kumo-install.XXXXXXXX")
trap 'rm -rf "$TMP"' EXIT HUP INT TERM
PREFIX="$TMP/prefix"
meson setup "$TMP/build-kumo" "$ROOT" --prefix="$PREFIX" --buildtype=debugoptimized
meson compile -C "$TMP/build-kumo"
meson install -C "$TMP/build-kumo"
LIBDIR=$(find "$PREFIX" -type f \( -name 'libkumo.so*' -o -name 'libkumo.dylib' -o -name 'libkumo.a' \) -print -quit)
if [ -z "$LIBDIR" ]; then echo 'ERROR: installed KUMO library not found' >&2; exit 1; fi
LIBDIR=$(dirname "$LIBDIR")
meson setup "$TMP/build-consumer" "$ROOT/tests/installed-consumer" -Dkumo_prefix="$PREFIX" --libdir="${LIBDIR#"$PREFIX"/}"
meson compile -C "$TMP/build-consumer"
if [ "$(uname -s)" = Darwin ]; then
  DYLD_LIBRARY_PATH="$LIBDIR${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" "$TMP/build-consumer/installed-kumo-consumer"
else
  LD_LIBRARY_PATH="$LIBDIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$TMP/build-consumer/installed-kumo-consumer"
fi

# Verify the installed CLI against the version declared by Meson.
if [ ! -x "$PREFIX/bin/kumo" ]; then
  echo 'ERROR: installed kumo CLI not found' >&2
  exit 1
fi

EXPECTED_VERSION=$(meson introspect "$TMP/build-kumo" --projectinfo |
  python3 -c 'import json,sys; print(json.load(sys.stdin)["version"])')

if [ "$(uname -s)" = Darwin ]; then
  ACTUAL_VERSION=$(DYLD_LIBRARY_PATH="$LIBDIR${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" "$PREFIX/bin/kumo" --version)
else
  ACTUAL_VERSION=$(LD_LIBRARY_PATH="$LIBDIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$PREFIX/bin/kumo" --version)
fi

case "$ACTUAL_VERSION" in
  *"$EXPECTED_VERSION"*) ;;
  *)
    echo "ERROR: installed CLI version mismatch" >&2
    echo "Expected: $EXPECTED_VERSION" >&2
    echo "Actual: $ACTUAL_VERSION" >&2
    exit 1
    ;;
esac

echo 'installed KUMO CLI OK'
