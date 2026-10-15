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

# Verify that the installed CLI is present and linked against the installed library.
if [ ! -x "$PREFIX/bin/kumo" ]; then
  echo 'ERROR: installed kumo CLI not found' >&2; exit 1
fi
if [ "$(uname -s)" = Darwin ]; then
  DYLD_LIBRARY_PATH="$LIBDIR${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" "$PREFIX/bin/kumo" --version | grep -q '1.3.0-dev'
else
  LD_LIBRARY_PATH="$LIBDIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$PREFIX/bin/kumo" --version | grep -q '1.3.0-dev'
fi
echo 'installed KUMO CLI OK'
