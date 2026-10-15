project('chem-kumo','cpp',version:'1.3.0-dev',meson_version:'>=0.63',default_options:['cpp_std=c++20','warning_level=2'])
# Single source of truth for version strings in the CLI and Cube headers.
add_project_arguments('-DCHEM_KUMO_VERSION="' + meson.project_version() + '"', language: 'cpp')
eigen_dep = dependency('eigen3', required: true)
zlib_dep = dependency('zlib', required: true)
# Optional OpenMP for independent grid blocks in density, derivative and MO evaluations.
# Fail explicitly if the user requests OpenMP but the toolchain lacks it.
openmp_enabled = get_option('openmp')
openmp_deps = []
if openmp_enabled
  openmp_deps += dependency('openmp', required: true)
  # Eigen must not start nested worker pools inside OpenMP grid workers.
  add_project_arguments('-DEIGEN_DONT_PARALLELIZE', language: 'cpp')
endif
message('OpenMP support requested: ' + (openmp_enabled ? 'yes' : 'no'))
benchmarks_enabled = get_option('benchmarks')
inc = include_directories('include')
kumo_lib = library('kumo', ['src/fchk.cpp','src/basis.cpp','src/basis_math.cpp','src/basis_evaluator.cpp','src/density.cpp','src/orbitals.cpp'],include_directories:inc,dependencies:[eigen_dep,zlib_dep] + openmp_deps,install:true)
kumo_dep = declare_dependency(link_with:kumo_lib,include_directories:inc,dependencies:[eigen_dep,zlib_dep] + openmp_deps)
meson.override_dependency('chem-kumo', kumo_dep)
# Installed consumers can discover the public headers and shared library.
# Eigen appears in public headers; zlib is an implementation dependency.
pkg = import('pkgconfig')
pkg.generate(
  kumo_lib,
  name: 'chem-kumo',
  filebase: 'chem-kumo',
  description: 'Ground-state electronic fields from Gaussian FCHK files',
  version: meson.project_version(),
  requires: ['eigen3'],
  requires_private: ['zlib'],
)
install_headers('include/chem/kumo.hpp',subdir:'chem')
install_headers('include/chem/kumo/density.hpp','include/chem/kumo/fchk.hpp','include/chem/kumo/basis.hpp','include/chem/kumo/basis_evaluator.hpp',subdir:'chem/kumo')
executable('kumo-example','examples/fchk_to_cube.cpp',dependencies:kumo_dep)
test('density',executable('kumo-density-test','tests/density_test.cpp',dependencies:kumo_dep),args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

cube_reference_test = executable('kumo-cube-reference-test', 'tests/cube_reference_test.cpp', dependencies:kumo_dep)
test('cube-reference-o2', cube_reference_test, args:[meson.current_source_dir()/'tests/data/o2.fchk.gz', meson.current_source_dir()/'tests/data/o2_density.cube'], timeout:300)
test('cube-reference-acrolein', cube_reference_test, args:[meson.current_source_dir()/'tests/data/acrolein.fchk.gz', meson.current_source_dir()/'tests/data/acrolein_density.cube'], timeout:300)

foreach component : ['alpha', 'beta', 'spin']
  test('cube-reference-o2-' + component, cube_reference_test,
    args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
          meson.current_source_dir()/('tests/data/o2_' + component + '.cube'), component], timeout:300)
endforeach

test('spin-densities', executable('kumo-spin-densities-test', 'tests/spin_densities_test.cpp', dependencies:kumo_dep), args:[meson.current_source_dir()/'tests/data/o2.fchk.gz', meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

test('selective-densities', executable('kumo-selective-densities-test',
  'tests/selective_densities_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])
if benchmarks_enabled
  executable('kumo-benchmark', 'benchmarks/density_benchmark.cpp', dependencies:kumo_dep)
endif

if benchmarks_enabled
  executable('kumo-profile', 'benchmarks/profile.cpp', dependencies:kumo_dep)
endif


test('ao-optimized', executable('kumo-ao-optimized-test', 'tests/ao_optimized_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz', meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])
if benchmarks_enabled
  executable('kumo-ao-benchmark', 'benchmarks/ao_benchmark.cpp', dependencies:kumo_dep)
endif

blocked_test = executable('kumo-blocked-test', 'tests/blocked_test.cpp', dependencies:kumo_dep)
test('blocked-densities', blocked_test,
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz', meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])
if benchmarks_enabled
  executable('kumo-block-benchmark', 'benchmarks/block_benchmark.cpp', dependencies:kumo_dep)
endif

if benchmarks_enabled
  executable('kumo-parallel-density-benchmark',
    'benchmarks/parallel_density_benchmark.cpp', dependencies:kumo_dep)
endif

if benchmarks_enabled
  executable('kumo-parallel-fields-benchmark',
    'benchmarks/parallel_fields_benchmark.cpp', dependencies:kumo_dep)
endif

# The default public path uses blocked Eigen contractions.  The original
# pointwise implementation is retained as an independent regression oracle.
test('default-blocked', executable('kumo-default-blocked-test',
  'tests/default_blocked_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

# Input-validation and public API invariants.
test('input-validation', executable('kumo-input-validation-test',
  'tests/input_validation_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz'])

# Public-header consumer smoke test, linked only through kumo_dep.
test('public-consumer', executable('kumo-public-consumer-test',
  'tests/consumer/main.cpp', dependencies:kumo_dep))

# Command-line interface, built entirely on the public library API.
kumo_cli = executable('kumo', 'src/kumo_main.cpp', dependencies:kumo_dep, install:true)

# Signed molecular-orbital amplitudes against Gaussian Cubegen references.
orbital_reference_test = executable('kumo-orbital-reference-test',
  'tests/orbital_reference_test.cpp', dependencies:kumo_dep)
foreach orbital : ['o2_alpha_mo_3', 'o2_alpha_mo_9', 'o2_alpha_mo_10',
                   'o2_beta_mo_3', 'o2_beta_mo_7', 'o2_beta_mo_8',
                   'acrolein_mo_13', 'acrolein_mo_15',
                   'acrolein_mo_16', 'acrolein_mo_18']
  parts = orbital.split('_')
  molecule = parts[0]
  channel = molecule == 'o2' ? parts[1] : 'alpha'
  index = parts[parts.length() - 1]
  test('orbital-reference-' + orbital, orbital_reference_test,
    args: [meson.current_source_dir() / ('tests/data/' + molecule + '.fchk.gz'),
           meson.current_source_dir() / ('tests/data/mo_reference_cubes/' + orbital + '.cube'),
           channel, index], timeout:300)
endforeach

kumo_cli_test = find_program('tests/cli_list_test.sh')
test('cli-list-restricted', kumo_cli_test, args:[kumo_cli, meson.current_source_dir()/'tests/data/acrolein.fchk.gz', 'restricted'])
test('cli-list-unrestricted', kumo_cli_test, args:[kumo_cli, meson.current_source_dir()/'tests/data/o2.fchk.gz', 'alpha'])

# Exercise the installed command-line interface on a real FCHK file.
kumo_cli_e2e = find_program('tests/cli_end_to_end.sh')
test('cli-end-to-end', kumo_cli_e2e,
  args:[kumo_cli, meson.current_source_dir()/'tests/data/o2.fchk.gz'], timeout:300)

# Spatial AO derivatives (analytic implementation versus numerical oracle).
test('ao-spatial-derivatives', executable('kumo-ao-spatial-derivatives-test',
  'tests/ao_derivatives_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

# Analytic spatial derivatives of signed molecular orbitals.
test('mo-spatial-derivatives', executable('kumo-mo-spatial-derivatives-test',
  'tests/mo_derivatives_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

# Analytic spatial derivatives of spin-resolved electron densities.
test('density-spatial-derivatives', executable('kumo-density-spatial-derivatives-test',
  'tests/density_derivatives_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

# Independent MO-occupation reconstruction and blocked-versus-pointwise derivatives.
test('density-mo-block-consistency', executable('kumo-density-mo-block-consistency-test',
  'tests/density_mo_consistency_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

# Optional profiling of spatial derivatives; not part of the CI test suite.
if benchmarks_enabled
  executable('kumo-spatial-benchmark', 'benchmarks/spatial_derivatives_benchmark.cpp', dependencies:kumo_dep)
endif

# Gzip regression tests use committed compressed references and materialize
# plain-text copies only inside a temporary test directory.
test('fchk-gzip', executable('kumo-fchk-gzip-test',
  'tests/fchk_gzip_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz',
        meson.current_source_dir()/'tests/data/hbc/hbc_b3lyp_631gdp.fchk.gz',
        meson.current_source_dir()/'tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz',
        meson.current_source_dir()/'tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz'])

# Compare density grids with 1/2/4 OpenMP workers (also runs in serial builds).
parallel_density_test = executable('kumo-parallel-density-test',
  'tests/parallel_density_test.cpp', dependencies:kumo_dep)
test('parallel-density', parallel_density_test,
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)
parallel_derivatives_test = executable('kumo-parallel-derivatives-test',
  'tests/parallel_derivatives_test.cpp', dependencies:kumo_dep)
test('parallel-derivatives', parallel_derivatives_test,
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# MO grid and spatial-derivative consistency across OpenMP worker counts.
parallel_orbitals_test = executable('kumo-parallel-orbitals-test',
  'tests/parallel_orbitals_test.cpp', dependencies:kumo_dep)
test('parallel-orbitals', parallel_orbitals_test,
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)
