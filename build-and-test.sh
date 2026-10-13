#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT}/build"
CLEAN=0

usage() {
  echo "Usage: $0 [--clean] [--help]"
  echo "  --clean  Remove build directory before configuring"
}

for arg in "$@"; do
  case "$arg" in
    --clean) CLEAN=1 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

for cmd in meson ninja; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: missing required command: $cmd" >&2
    exit 1
  fi
done

if [[ ! -f "${ROOT}/meson.build" ]]; then
  echo "ERROR: meson.build not found in ${ROOT}" >&2
  echo "Place this script in the chem-kumo project root." >&2
  exit 1
fi

if (( CLEAN )); then
  echo "==> Removing ${BUILD_DIR}"
  rm -rf -- "$BUILD_DIR"
fi

if [[ -f "${BUILD_DIR}/build.ninja" ]]; then
  echo "==> Reconfiguring Meson"
  meson setup --reconfigure "$BUILD_DIR" "$ROOT"
else
  echo "==> Configuring Meson"
  meson setup "$BUILD_DIR" "$ROOT" --buildtype=debugoptimized
fi

echo "==> Compiling chem-kumo"
meson compile -C "$BUILD_DIR"

echo "==> Running all tests"
meson test -C "$BUILD_DIR" --print-errorlogs --verbose

echo "==> SUCCESS: build and tests completed"
