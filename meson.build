project('chem-kumo','cpp',version:'1.1.0',meson_version:'>=0.63',default_options:['cpp_std=c++20','warning_level=2'])
eigen_dep = dependency('eigen3', required: true)
benchmarks_enabled = get_option('benchmarks')
inc = include_directories('include')
kumo_lib = library('kumo', ['src/fchk.cpp','src/basis.cpp','src/basis_math.cpp','src/basis_evaluator.cpp','src/density.cpp','src/orbitals.cpp'],include_directories:inc,dependencies:eigen_dep,install:true)
kumo_dep = declare_dependency(link_with:kumo_lib,include_directories:inc,dependencies:eigen_dep)
meson.override_dependency('chem-kumo', kumo_dep)
install_headers('include/chem/kumo.hpp',subdir:'chem')
install_headers('include/chem/kumo/density.hpp','include/chem/kumo/fchk.hpp','include/chem/kumo/basis.hpp','include/chem/kumo/basis_evaluator.hpp',subdir:'chem/kumo')
executable('kumo-example','examples/fchk_to_cube.cpp',dependencies:kumo_dep)
test('density',executable('kumo-density-test','tests/density_test.cpp',dependencies:kumo_dep),args:[meson.current_source_dir()/'tests/data/o2.fchk',meson.current_source_dir()/'tests/data/acrolein.fchk'])

cube_reference_test = executable('kumo-cube-reference-test', 'tests/cube_reference_test.cpp', dependencies:kumo_dep)
test('cube-reference-o2', cube_reference_test, args:[meson.current_source_dir()/'tests/data/o2.fchk', meson.current_source_dir()/'tests/data/o2_density.cube'], timeout:300)
test('cube-reference-acrolein', cube_reference_test, args:[meson.current_source_dir()/'tests/data/acrolein.fchk', meson.current_source_dir()/'tests/data/acrolein_density.cube'], timeout:300)

foreach component : ['alpha', 'beta', 'spin']
  test('cube-reference-o2-' + component, cube_reference_test,
    args:[meson.current_source_dir()/'tests/data/o2.fchk',
          meson.current_source_dir()/('tests/data/o2_' + component + '.cube'), component], timeout:300)
endforeach

test('spin-densities', executable('kumo-spin-densities-test', 'tests/spin_densities_test.cpp', dependencies:kumo_dep), args:[meson.current_source_dir()/'tests/data/o2.fchk', meson.current_source_dir()/'tests/data/acrolein.fchk'])

test('selective-densities', executable('kumo-selective-densities-test',
  'tests/selective_densities_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk',
        meson.current_source_dir()/'tests/data/acrolein.fchk'])
if benchmarks_enabled
  executable('kumo-benchmark', 'benchmarks/density_benchmark.cpp', dependencies:kumo_dep)
endif

if benchmarks_enabled
  executable('kumo-profile', 'benchmarks/profile.cpp', dependencies:kumo_dep)
endif


test('ao-optimized', executable('kumo-ao-optimized-test', 'tests/ao_optimized_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk', meson.current_source_dir()/'tests/data/acrolein.fchk'])
if benchmarks_enabled
  executable('kumo-ao-benchmark', 'benchmarks/ao_benchmark.cpp', dependencies:kumo_dep)
endif

blocked_test = executable('kumo-blocked-test', 'tests/blocked_test.cpp', dependencies:kumo_dep)
test('blocked-densities', blocked_test,
  args:[meson.current_source_dir()/'tests/data/o2.fchk', meson.current_source_dir()/'tests/data/acrolein.fchk'])
if benchmarks_enabled
  executable('kumo-block-benchmark', 'benchmarks/block_benchmark.cpp', dependencies:kumo_dep)
endif

# The default public path uses blocked Eigen contractions.  The original
# pointwise implementation is retained as an independent regression oracle.
test('default-blocked', executable('kumo-default-blocked-test',
  'tests/default_blocked_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk',
        meson.current_source_dir()/'tests/data/acrolein.fchk'])

# Input-validation and public API invariants.
test('input-validation', executable('kumo-input-validation-test',
  'tests/input_validation_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk'])

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
    args: [meson.current_source_dir() / ('tests/data/' + molecule + '.fchk'),
           meson.current_source_dir() / ('tests/data/mo_reference_cubes/' + orbital + '.cube'),
           channel, index], timeout:300)
endforeach

kumo_cli_test = find_program('tests/cli_list_test.sh')
test('cli-list-restricted', kumo_cli_test, args:[kumo_cli, meson.current_source_dir()/'tests/data/acrolein.fchk', 'restricted'])
test('cli-list-unrestricted', kumo_cli_test, args:[kumo_cli, meson.current_source_dir()/'tests/data/o2.fchk', 'alpha'])

# Exercise the installed command-line interface on a real FCHK file.
kumo_cli_e2e = find_program('tests/cli_end_to_end.sh')
test('cli-end-to-end', kumo_cli_e2e,
  args:[kumo_cli, meson.current_source_dir()/'tests/data/o2.fchk'], timeout:300)
