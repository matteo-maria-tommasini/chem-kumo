#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CLEAN=0
OPENMP=0
BUILD_TYPE=debugoptimized

usage() {
  cat <<'HELP'
Usage: ./build-and-test.sh [--clean] [--openmp] [--release] [--help]
  --clean    Recreate only the selected build directory
  --openmp   Enable OpenMP (MacPorts Clang 22 on macOS, default GCC on Linux)
  --release  Use Meson release optimization (default: debugoptimized)

Build directories: build/ (serial), build-openmp/ (OpenMP).
HELP
}

for arg in "$@"; do
  case "$arg" in
    --clean) CLEAN=1 ;;
    --openmp) OPENMP=1 ;;
    --release) BUILD_TYPE=release ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

for cmd in meson ninja; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: missing $cmd" >&2; exit 1; }
done
[[ -f "$ROOT/meson.build" ]] || { echo "ERROR: meson.build missing in $ROOT" >&2; exit 1; }

BUILD_DIR="$ROOT/build"
OPENMP_SETTING=false
MESON_NATIVE=()
if (( OPENMP )); then
  BUILD_DIR="$ROOT/build-openmp"
  OPENMP_SETTING=true
  if [[ "$(uname -s)" == Darwin ]]; then
    NATIVE_FILE="$ROOT/meson/native/macports-clang22.ini"
    [[ -f "$NATIVE_FILE" ]] || { echo "ERROR: missing $NATIVE_FILE" >&2; exit 1; }
    for compiler in /opt/local/bin/clang-mp-22 /opt/local/bin/clang++-mp-22; do
      [[ -x "$compiler" ]] || { echo "ERROR: missing compiler $compiler" >&2; exit 1; }
    done
    MESON_NATIVE=(--native-file "$NATIVE_FILE")
  fi
fi

if (( CLEAN )); then
  echo "==> Removing $BUILD_DIR"
  rm -rf -- "$BUILD_DIR"
fi

if [[ -f "$BUILD_DIR/build.ninja" ]]; then
  echo "==> Reconfiguring $BUILD_DIR (OpenMP=$OPENMP_SETTING, buildtype=$BUILD_TYPE)"
  meson setup --reconfigure "$BUILD_DIR" "$ROOT" \
    "-Dopenmp=$OPENMP_SETTING" "-Dbuildtype=$BUILD_TYPE"
else
  if [[ -e "$BUILD_DIR" ]]; then
    echo "ERROR: $BUILD_DIR exists but is not a configured Ninja build; remove it with --clean" >&2
    exit 1
  fi
  echo "==> Configuring $BUILD_DIR (OpenMP=$OPENMP_SETTING, buildtype=$BUILD_TYPE)"
  meson setup "$BUILD_DIR" "$ROOT" "${MESON_NATIVE[@]}" \
    "-Dopenmp=$OPENMP_SETTING" "-Dbuildtype=$BUILD_TYPE"
fi

echo "==> Compiling chem-kumo"
meson compile -C "$BUILD_DIR"

echo "==> Running all tests"
CI=1 meson test -C "$BUILD_DIR" --print-errorlogs

echo "==> SUCCESS: build and tests completed"
