# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project

chem-kumo is a C++20 library (`libkumo`) and CLI (`kumo`) that reconstructs
ground-state SCF electron densities, density derivatives and molecular
orbitals from Gaussian formatted checkpoint files (`.fchk`, `.fchk.gz`) and
writes/reads scalar Gaussian Cube files (also gzip-compressed).
Dependencies: Eigen 3.4 (public, appears in headers), zlib (private),
optional OpenMP. Public GitHub repo: `matteo-maria-tommasini/chem-kumo`.

## Layout

- `include/chem/kumo.hpp`, `include/chem/kumo/*.hpp` — public API (installed).
- `src/*.cpp` — library implementation; `src/*.hpp` are internal headers
  (not installed). `src/kumo_main.cpp` is the CLI.
- `tests/` — one executable per `*_test.cpp`, plus shell tests
  (`cli_*.sh`, `check-installed-consumer.sh`). Reference data in
  `tests/data/` (Gaussian `.fchk.gz` and reference Cube files; provenance in
  `tests/REFERENCE_CUBES.md`).
- `benchmarks/` — optional, built only with `-Dbenchmarks=true`.
- `docs/API.md`, `docs/OPENMP.md`, `docs/mathematics/*.tex` — documentation.
- `meson.build` — the single build file; every test is registered there.

## Build and test

```sh
./build-and-test.sh --clean              # serial build in build/
./build-and-test.sh --clean --openmp     # OpenMP build in build-openmp/
meson compile -C build && meson test -C build --print-errorlogs
./tests/check-installed-consumer.sh      # install + external consumer check
```

On macOS, Eigen/Meson/Ninja come from MacPorts (`/opt/local`); OpenMP builds
use `meson/native/macports-clang22.ini`. Build dirs (`build/`, `build-*/`)
and generated `*.cube` outside `tests/data/` are git-ignored.

Always run the full test suite before proposing a commit. Reference tests
compare against Gaussian output, so do not loosen numerical tolerances to
make a test pass without explicit approval.

## Conventions

- The version lives only in `project(... version: ...)` in `meson.build`; it
  is passed to code as `CHEM_KUMO_VERSION`. Tests derive the expected version
  from Meson — never hard-code it. On a release, also update the README title
  and `RELEASE_NOTES.md`.
- New public headers must be added to `install_headers(...)` in
  `meson.build`; new sources to the `kumo` library source list.
- New tests: add `tests/<name>_test.cpp` and register a `test(...)` in
  `meson.build` (with a `timeout` for slow reference tests).
- Validate inputs strictly and fail with clear exceptions (fail-closed);
  keep serial and OpenMP results numerically identical in point ordering.
- Compiler warnings are enabled (`warning_level=2`); keep the build warning-free.

## Git and GitHub workflow

- `main` is the default branch; develop on `feature/<topic>` branches and
  merge via pull request. CI (`.github/workflows/ci.yml`) builds and tests on
  ubuntu-latest and macos-latest and runs the installed-consumer check.
- Commit messages: Conventional Commits style (`feat:`, `fix:`, `test:`,
  `docs:`, `build:`, `ci:`, `refactor:`), imperative, concise.
- Releases are annotated tags `vX.Y.Z` on `main`.
- Use `git` and `gh` from the shell. Ask before pushing, opening/merging PRs,
  tagging or creating releases — the repository is public.
