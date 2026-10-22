project('chem-kumo','cpp',version:'1.3.0',meson_version:'>=0.63',default_options:['cpp_std=c++20','warning_level=2'])
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
kumo_lib = library('kumo', ['src/fchk.cpp','src/basis.cpp','src/basis_math.cpp','src/basis_evaluator.cpp','src/density.cpp','src/orbitals.cpp','src/evaluator.cpp'],include_directories:inc,dependencies:[eigen_dep,zlib_dep] + openmp_deps,install:true)
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
install_headers('include/chem/kumo/density.hpp','include/chem/kumo/fchk.hpp','include/chem/kumo/basis.hpp','include/chem/kumo/basis_evaluator.hpp','include/chem/kumo/evaluator.hpp',subdir:'chem/kumo')
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

if benchmarks_enabled
  executable('kumo-backend-comparison',
    'benchmarks/backend_comparison_benchmark.cpp', dependencies:kumo_dep)
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

# Stage 01C: AO-P versus AO-MO derivative benchmark.
if benchmarks_enabled
  executable(
    'kumo-derivative-backend-comparison',
    'benchmarks/derivative_backend_benchmark.cpp',
    dependencies: kumo_dep
  )
endif

test('evaluator-api', executable('kumo-evaluator-api-test', 'tests/evaluator_api_test.cpp', dependencies:kumo_dep), args:[meson.current_source_dir()/'tests/data/o2.fchk.gz'])

# Stage 02B.1: fused/individual evaluator requests against validated 1.3 kernels.
test('evaluator-fusion', executable('kumo-evaluator-fusion-test',
  'tests/evaluator_fusion_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 02C.1: internal dependency planning, without changing numerical kernels.
test('evaluation-plan', executable('kumo-evaluation-plan-test',
  ['tests/evaluation_plan_test.cpp'], include_directories: include_directories('src'),
  dependencies: kumo_dep))

# Stage 02C.2c: direct numerical regression against the legacy AO-P derivative kernel.
test('shared-ao-validation', executable('kumo-shared-ao-validation-test',
  'tests/shared_ao_validation_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Optional single-thread algorithmic comparison; run manually, not in CI.
if benchmarks_enabled
  executable('kumo-shared-ao-benchmark', 'benchmarks/shared_ao_benchmark.cpp', dependencies:kumo_dep)
endif

# Stage 02C.3: first-order AO path against full second-order contractions.
test('selective-ao-derivatives', executable('kumo-selective-ao-derivatives-test',
  'tests/selective_ao_derivatives_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 02C.3b: benchmark selective first-order versus second-order AO derivatives.
if benchmarks_enabled
  executable('kumo-selective-ao-benchmark',
    'benchmarks/selective_ao_benchmark.cpp',
    dependencies: kumo_dep)
endif

# Stage 02D.1: validate AO-MO eligibility without changing production contractions.
test('backend-selection', executable('kumo-backend-selection-test',
  'tests/backend_selection_test.cpp', include_directories: include_directories('src'),
  dependencies: kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'])

# Stage 02D.2a: experimental AO-MO kernel, independent of Evaluator selection.
test('ao-mo-kernel', executable('kumo-ao-mo-kernel-test',
  'tests/ao_mo_kernel_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 02D.2a: manual three-way AO-MO numerical diagnostic.
if benchmarks_enabled
  executable('kumo-ao-mo-diagnostic',
    'tests/ao_mo_diagnostic.cpp', dependencies: kumo_dep)
endif

# Stage 02D.2b: explicit AO_P / AO_MO Evaluator backend selection.
test('evaluator-backend', executable('kumo-evaluator-backend-test',
  'tests/evaluator_backend_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 02D.2c: manual end-to-end AO_P vs AO_MO comparison through Evaluator.
if benchmarks_enabled
  executable('kumo-evaluator-backend-benchmark',
    'benchmarks/evaluator_backend_benchmark.cpp', dependencies:kumo_dep)
endif

# Stage 02D.3a: deterministic, conservative automatic backend selection.
test('automatic-backend', executable('kumo-automatic-backend-test',
  'tests/automatic_backend_test.cpp', include_directories: include_directories('src'),
  dependencies: kumo_dep))

# Stage 02D.3b: measure AUTO decisions against both explicit backends.
if benchmarks_enabled
  executable('kumo-automatic-backend-benchmark',
    'benchmarks/automatic_backend_benchmark.cpp',
    include_directories: include_directories('src'),
    dependencies: kumo_dep)
endif

# Stage 03A: analytic AO overlaps and Gaussian MO metric regression.
test('analytic-overlap', executable('kumo-analytic-overlap-test',
  'tests/analytic_overlap_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 03A.2: independent Gauss-Hermite primitive-product overlap reference.
test('analytic-overlap-reference', executable('kumo-analytic-overlap-reference-test',
  'tests/analytic_overlap_reference_test.cpp', dependencies: kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 03B.1: expose the stored spin-resolved FCHK 1-RDM without duplicating state.
test('one-particle-density', executable('kumo-one-particle-density-test',
  'tests/one_particle_density_test.cpp', dependencies: kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 03B.2: S-metric subspace projectors (factorized representation).
test('subspace-projector', executable('kumo-subspace-projector-test',
  'tests/subspace_projector_test.cpp', dependencies: kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 03B.2a: metric projector stress tests.
test('subspace-projector-stress',
  executable('kumo-subspace-projector-stress-test',
    'tests/subspace_projector_stress_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.1: quadrature-weighted AO overlap accumulation (no grid policy).
test('weighted-overlap', executable('kumo-weighted-overlap-test',
  'tests/weighted_overlap_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 03C.2: diagnostic convergence of spatial quadrature toward analytic AO overlap.
test('weighted-overlap-quadrature', executable('kumo-weighted-overlap-quadrature-test',
  'tests/weighted_overlap_quadrature_test.cpp', dependencies:kumo_dep),
  args:[meson.current_source_dir()/'tests/data/o2.fchk.gz',
        meson.current_source_dir()/'tests/data/acrolein.fchk.gz'], timeout:300)

# Stage 03C.3: independent atomic quadrature validation.
test('atomic-radial-angular-quadrature',
  executable('kumo-atomic-radial-angular-quadrature-test',
    'tests/atomic_radial_angular_quadrature_test.cpp',
    dependencies: kumo_dep),
  timeout:300)

# Stage 03C.4: higher-order single-centre quadrature convergence.
test('atomic-radial-angular-convergence', executable(
  'kumo-atomic-radial-angular-convergence-test',
  'tests/atomic_radial_angular_convergence_test.cpp',
  dependencies: kumo_dep), timeout:300)

# Stage 03C.5: end-to-end single-center AO overlap quadrature.
test('weighted-overlap-single-center',
  executable('kumo-weighted-overlap-single-center-test',
    'tests/weighted_overlap_single_center_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.6: Becke multicenter overlap diagnostic.
test('weighted-overlap-multicenter', executable(
  'kumo-weighted-overlap-multicenter-test',
  'tests/weighted_overlap_multicenter_test.cpp',
  dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.8a: analytic local four-Gaussian primitive integrals.
test('local-four-gaussian-integrals', executable(
  'kumo-local-four-gaussian-integrals-test',
  'tests/local_four_gaussian_integral_test.cpp',
  dependencies: kumo_dep), timeout: 300)
# Stage 03C.8b: contracted Cartesian/pure local shell quartets.
test('local-shell-quartet', executable('kumo-local-shell-quartet-test',
  'tests/local_shell_quartet_test.cpp', dependencies: kumo_dep),
  args: [meson.current_source_dir()/'tests/data/o2.fchk.gz',
         meson.current_source_dir()/'tests/data/acrolein.fchk.gz'],
  timeout: 300)

# Stage 03C.8b.1: independent contracted shell quartet reference.
test('local-shell-quartet-reference',
  executable('kumo-local-shell-quartet-reference-test',
    'tests/local_shell_quartet_reference_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.8c: analytic density-weighted AO overlap.
test('local-density-overlap',
  executable('kumo-local-density-overlap-test',
    'tests/local_density_overlap_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.8d: independent numerical density-overlap reference.
test('local-density-overlap-reference',
  executable('kumo-local-density-overlap-reference-test',
    'tests/local_density_overlap_reference_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.9a: symmetry-reduced analytic density overlap.
test('local-density-overlap-symmetric',
  executable('kumo-local-density-overlap-symmetric-test',
    'tests/local_density_overlap_symmetric_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.9b: Gaussian primitive-pair profiling.
executable('kumo-local-gaussian-pair-profile',
  'tests/local_gaussian_pair_profile.cpp',
  dependencies: kumo_dep)

# Stage 03C.9d: Gaussian-pair shell quartet validation.
test('local-shell-pair-integrals',
  executable('kumo-local-shell-pair-integrals-test',
    'tests/local_shell_pair_integrals_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.9e: symmetric density overlap with Gaussian-pair intermediates.
test('local-density-overlap-symmetric-paired',
  executable('kumo-local-density-overlap-symmetric-paired-test',
    'tests/local_density_overlap_symmetric_paired_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.10a: analytic AO product norms and Cauchy-Schwarz bounds.
test('local-pair-bounds',
  executable('kumo-local-pair-bounds-test',
    'tests/local_pair_bounds_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.10b: density-aware rigorous screening diagnostic.
test('local-screening-diagnostic',
  executable('kumo-local-screening-diagnostic-test',
    'tests/local_screening_diagnostic_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.10f: contracted shell-pair Gaussian envelopes.
test('local-shell-pair-envelope',
  executable('kumo-local-shell-pair-envelope-test',
    'tests/local_shell_pair_envelope_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.10g: comparison of shell-level Gaussian bounds.
test('local-shell-bound-comparison',
  executable('kumo-local-shell-bound-comparison-test',
    'tests/local_shell_bound_comparison_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.11a: error-budgeted screened density overlap.
test('local-density-overlap-screened',
  executable('kumo-local-density-overlap-screened-test',
    'tests/local_density_overlap_screened_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.11b: reusable density-overlap preparation.
test('local-density-overlap-plan',
  executable('kumo-local-density-overlap-plan-test',
    'tests/local_density_overlap_plan_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.11c: memory and scalability diagnostics.
executable('kumo-local-density-overlap-memory-test',
  'tests/local_density_overlap_memory_test.cpp',
  dependencies: kumo_dep)

# Stage 03C.11d: streaming CS screening diagnostic.
executable('kumo-local-screening-streaming-test',
  'tests/local_screening_streaming_test.cpp',
  dependencies: kumo_dep)

# Stage 03C.11e: bounded-memory streaming density overlap.
test('local-density-overlap-streaming',
  executable('kumo-local-density-overlap-streaming-test',
    'tests/local_density_overlap_streaming_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.12a: hierarchical screening diagnostic.
executable('kumo-local-screening-hierarchical-test',
  'tests/local_screening_hierarchical_test.cpp',
  dependencies: kumo_dep)

# Stage 03C.12b: lazy hierarchical screening diagnostic.
executable('kumo-local-screening-lazy-test',
  'tests/local_screening_lazy_test.cpp',
  dependencies: kumo_dep)

# Stage 03C.12c: importance-ordered screening diagnostic.
executable('kumo-local-screening-importance-test',
  'tests/local_screening_importance_test.cpp',
  dependencies: kumo_dep)

# Stage 03C.12d: importance-ordered streaming contraction.
test('local-density-overlap-importance',
  executable('kumo-local-density-overlap-importance-test',
    'tests/local_density_overlap_importance_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.13a: independent AO symmetry regression.
symmetry_regression_exe = executable(
  'kumo-local-density-overlap-symmetry-regression-test',
  'tests/local_density_overlap_symmetry_regression_test.cpp',
  dependencies: kumo_dep
)

test(
  'local-density-overlap-symmetry-regression',
  symmetry_regression_exe,
  args: [
    meson.current_source_dir() / 'tests/data/o2.fchk.gz'
  ],
  timeout: 300
)

# Stage 03C.13c: shared-quartet multichannel contraction.
test('local-density-overlap-multichannel',
  executable('kumo-local-density-overlap-multichannel-test',
    'tests/local_density_overlap_multichannel_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.13c1: multichannel edge-case regression.
test('local-density-overlap-multichannel-edge',
  executable('kumo-local-density-overlap-multichannel-edge-test',
    'tests/local_density_overlap_multichannel_edge_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.13d: shared-quartet screening with independent budgets.
test('local-density-overlap-multichannel-screened',
  executable('kumo-local-density-overlap-multichannel-screened-test',
    'tests/local_density_overlap_multichannel_screened_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.13e1: multichannel importance-order diagnostic.
executable('kumo-local-screening-multichannel-importance-test',
  'tests/local_screening_multichannel_importance_test.cpp',
  dependencies: kumo_dep)

# Stage 03C.13e2: importance-ordered multichannel contraction.
test('local-density-overlap-multichannel-importance',
  executable('kumo-local-density-overlap-multichannel-importance-test',
    'tests/local_density_overlap_multichannel_importance_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.14a: large-system screening and AO-norm profile.
executable('kumo-local-screening-large-system-profile',
  'tests/local_screening_large_system_profile.cpp',
  dependencies: kumo_dep)

# Stage 03C.14c: AO product norm kernel profiling.
executable('kumo-local-ao-norm-profile',
  'tests/local_ao_norm_profile.cpp',
  dependencies: kumo_dep)

# Stage 03C.14d1-14e: AO norm regression and primitive symmetry.
test('local-ao-norm-regression',
  executable('kumo-local-ao-norm-regression-test',
    'tests/local_ao_norm_regression_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz',
    meson.current_source_dir()/'tests/data/acrolein.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.15a: directed aggregate-tail screening diagnostic.
executable('kumo-local-screening-block-tail-diagnostic',
  'tests/local_screening_block_tail_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.15b1: bidirectional cutoff diagnostic.
executable('kumo-local-screening-bidirectional-cutoff-diagnostic',
  'tests/local_screening_bidirectional_cutoff_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.15b2: hierarchical multichannel contraction.
test('local-density-overlap-multichannel-hierarchical',
  executable('kumo-local-density-overlap-multichannel-hierarchical-test',
    'tests/local_density_overlap_multichannel_hierarchical_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.15b3: hierarchical edge regression and benchmark.
test('local-density-overlap-hierarchical-benchmark',
  executable('kumo-local-density-overlap-hierarchical-benchmark-test',
    'tests/local_density_overlap_hierarchical_benchmark_test.cpp',
    dependencies: kumo_dep),
  args: [
    meson.current_source_dir()/'tests/data/o2.fchk.gz'
  ],
  timeout: 300)

# Stage 03C.16a1: spatial Gaussian envelope aggregation.
executable('kumo-local-spatial-envelope-aggregation-diagnostic',
  'tests/local_spatial_envelope_aggregation_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.16a3: guarded spatial aggregate bound.
test('local-spatial-aggregate-guard',
  executable('kumo-local-spatial-aggregate-guard-test',
    'tests/local_spatial_aggregate_guard_test.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16b1: adaptive spatial envelope hierarchy.
test('local-spatial-adaptive',
  executable('kumo-local-spatial-adaptive-diagnostic',
    'tests/local_spatial_adaptive_diagnostic.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16b2: persistent spatial envelope hierarchy.
test('local-spatial-persistent',
  executable('kumo-local-spatial-persistent-diagnostic',
    'tests/local_spatial_persistent_diagnostic.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16b3: recursive/persistent tree equivalence.
test('local-spatial-tree-equivalence',
  executable('kumo-local-spatial-tree-equivalence-test',
    'tests/local_spatial_tree_equivalence_test.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16b4: persistent spatial tree reuse benchmark.
executable('kumo-local-spatial-tree-reuse-benchmark',
  'tests/local_spatial_tree_reuse_benchmark.cpp',
  dependencies: kumo_dep)

# Stage 03C.16b5: persistent spatial tree depth/budget sweep.
executable('kumo-local-spatial-depth-budget-sweep',
  'tests/local_spatial_depth_budget_sweep.cpp',
  dependencies: kumo_dep)

# Stage 03C.16b6: spatial decision and budget safety.
test('local-spatial-decision-safety',
  executable('kumo-local-spatial-decision-safety-test',
    'tests/local_spatial_decision_safety_test.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16b7: extended spatial decision regression.
test('local-spatial-decision-extended',
  executable('kumo-local-spatial-decision-extended-test',
    'tests/local_spatial_decision_extended_test.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16c1: hybrid CS/spatial bound diagnostic.
executable('kumo-local-hybrid-cs-spatial-diagnostic',
  'tests/local_hybrid_cs_spatial_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.16c2: distance-stratified hybrid bound diagnostic.
executable('kumo-local-hybrid-cs-spatial-distance-diagnostic',
  'tests/local_hybrid_cs_spatial_distance_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.16c3: physical-density-weighted hybrid diagnostic.
executable('kumo-local-hybrid-density-weighted-diagnostic',
  'tests/local_hybrid_density_weighted_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.16c4: sampled hybrid cumulative-budget regression.
test('local-hybrid-budget-regression',
  executable('kumo-local-hybrid-budget-regression-test',
    'tests/local_hybrid_budget_regression_test.cpp',
    dependencies: kumo_dep),
  args: [meson.current_source_dir()/'tests/data/o2.fchk.gz'],
  timeout: 300)

# Stage 03C.16c5: CS vs hybrid cumulative-budget comparison.
test('local-hybrid-cs-budget-comparison',
  executable('kumo-local-hybrid-cs-budget-comparison-test',
    'tests/local_hybrid_cs_budget_comparison_test.cpp',
    dependencies: kumo_dep),
  args: [meson.current_source_dir()/'tests/data/o2.fchk.gz'],
  timeout: 300)

# Stage 03C.16c6: CS/hybrid screening cost diagnostic.
executable('kumo-local-hybrid-cost-diagnostic',
  'tests/local_hybrid_cost_diagnostic.cpp',
  dependencies: kumo_dep)

# Stage 03C.16c7: local budget accounting.
test('local-sparse-budget',
  executable('kumo-local-sparse-budget-test',
    'tests/local_sparse_budget_test.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16c7b: budget rejection and transaction atomicity.
test('local-sparse-budget-rejection',
  executable('kumo-local-sparse-budget-rejection-test',
    'tests/local_sparse_budget_rejection_test.cpp',
    dependencies: kumo_dep),
  timeout: 120)

# Stage 03C.16c7c: dense/sparse molecular budget equivalence and timings.
executable('kumo-local-hybrid-dense-sparse-equivalence',
  'tests/local_hybrid_dense_sparse_equivalence.cpp',
  dependencies: kumo_dep)

# Stage 03C.16c8: independent CS/hybrid contraction benchmark.
executable('kumo-local-hybrid-real-contraction-benchmark',
  'tests/local_hybrid_real_contraction_benchmark.cpp',
  dependencies: kumo_dep)

# Stage 03C.16c9: repeated CS/hybrid contraction benchmark.
executable('kumo-local-hybrid-repeated-benchmark',
  'tests/local_hybrid_repeated_benchmark.cpp',
  dependencies: kumo_dep)
