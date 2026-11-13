project(
    'chem-kumo',
    'cpp',
    version: '1.9.0',
    meson_version: '>=0.63',
    default_options: ['cpp_std=c++20', 'warning_level=2'],
)
# Single source of truth for version strings in the CLI and Cube headers.
add_project_arguments('-DCHEM_KUMO_VERSION="' + meson.project_version() + '"', language: 'cpp')
# Eigen is included as a system dependency (-isystem), so that compiler
# warnings raised inside Eigen's own headers, such as GCC's known
# -Wmaybe-uninitialized false positives in Eigen 3.4, are not reported.
# Warnings in this project's code are unaffected.
eigen_dep = dependency('eigen3', required: true, include_type: 'system')
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
kumo_lib = library(
    'kumo',
    [
        'src/fchk.cpp',
        'src/basis.cpp',
        'src/basis_math.cpp',
        'src/basis_evaluator.cpp',
        'src/boys_function.cpp',
        'src/electrostatic_potential.cpp',
        'src/density.cpp',
        'src/cube.cpp',
        'src/orbitals.cpp',
        'src/evaluator.cpp',
        'src/gzip_input.cpp',
    ],
    include_directories: inc,
    dependencies: [eigen_dep, zlib_dep] + openmp_deps,
    install: true,
)
kumo_dep = declare_dependency(
    link_with: kumo_lib,
    include_directories: inc,
    dependencies: [eigen_dep, zlib_dep] + openmp_deps,
)
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
install_headers('include/chem/kumo.hpp', subdir: 'chem')
install_headers(
    'include/chem/kumo/density.hpp',
    'include/chem/kumo/fchk.hpp',
    'include/chem/kumo/basis.hpp',
    'include/chem/kumo/basis_evaluator.hpp',
    'include/chem/kumo/evaluator.hpp',
    'include/chem/kumo/field_sampler.hpp',
    'include/chem/kumo/field_integrals.hpp',
    'include/chem/kumo/cube.hpp',
    'include/chem/kumo/subspace_projector.hpp',
    'include/chem/kumo/weighted_overlap.hpp',
    subdir: 'chem/kumo',
)
executable('kumo-example', 'examples/fchk_to_cube.cpp', dependencies: kumo_dep)
test(
    'density',
    executable('kumo-density-test', 'tests/density_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

cube_reference_test = executable(
    'kumo-cube-reference-test',
    'tests/cube_reference_test.cpp',
    dependencies: kumo_dep,
)
test(
    'cube-reference-o2',
    cube_reference_test,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/o2_density.cube',
    ],
    timeout: 300,
)
test(
    'cube-reference-acrolein',
    cube_reference_test,
    args: [
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein_density.cube',
    ],
    timeout: 300,
)

foreach component : ['alpha', 'beta', 'spin']
    test(
        'cube-reference-o2-' + component,
        cube_reference_test,
        args: [meson.current_source_dir() / 'tests/data/o2.fchk.gz', meson.current_source_dir() / (
            'tests/data/o2_' + component + '.cube'
        ), component],
        timeout: 300,
    )
endforeach

test(
    'spin-densities',
    executable('kumo-spin-densities-test', 'tests/spin_densities_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

test(
    'selective-densities',
    executable(
        'kumo-selective-densities-test',
        'tests/selective_densities_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

density_screening_test_exe = executable(
    'kumo-density-screening-test',
    'tests/density_screening_test.cpp',
    dependencies: kumo_dep,
)
test(
    'density-screening',
    density_screening_test_exe,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 120,
)
# The large HBC case is in the 'slow' suite, which the sanitizer CI job skips.
test(
    'density-screening-hbc',
    density_screening_test_exe,
    args: [meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz'],
    suite: 'slow',
    timeout: 120,
)
if benchmarks_enabled
    executable('kumo-benchmark', 'benchmarks/density_benchmark.cpp', dependencies: kumo_dep)
endif

if benchmarks_enabled
    executable('kumo-profile', 'benchmarks/profile.cpp', dependencies: kumo_dep)
endif


test(
    'ao-optimized',
    executable('kumo-ao-optimized-test', 'tests/ao_optimized_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)
if benchmarks_enabled
    executable('kumo-ao-benchmark', 'benchmarks/ao_benchmark.cpp', dependencies: kumo_dep)
endif

blocked_test = executable('kumo-blocked-test', 'tests/blocked_test.cpp', dependencies: kumo_dep)
test(
    'blocked-densities',
    blocked_test,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)
if benchmarks_enabled
    executable('kumo-block-benchmark', 'benchmarks/block_benchmark.cpp', dependencies: kumo_dep)
endif

if benchmarks_enabled
    executable(
        'kumo-parallel-density-benchmark',
        'benchmarks/parallel_density_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

if benchmarks_enabled
    executable(
        'kumo-parallel-fields-benchmark',
        'benchmarks/parallel_fields_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

if benchmarks_enabled
    executable(
        'kumo-backend-comparison',
        'benchmarks/backend_comparison_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

# The default public path uses blocked Eigen contractions.  The original
# pointwise implementation is retained as an independent regression oracle.
test(
    'default-blocked',
    executable('kumo-default-blocked-test', 'tests/default_blocked_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Input-validation and public API invariants.
test(
    'input-validation',
    executable(
        'kumo-input-validation-test',
        'tests/input_validation_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [meson.current_source_dir() / 'tests/data/o2.fchk.gz'],
)

# Public-header consumer smoke test, linked only through kumo_dep.
test(
    'public-consumer',
    executable('kumo-public-consumer-test', 'tests/consumer/main.cpp', dependencies: kumo_dep),
)

# Command-line interface, built entirely on the public library API.
kumo_cli = executable('kumo', 'src/kumo_main.cpp', dependencies: kumo_dep, install: true)

# Signed molecular-orbital amplitudes against Gaussian Cubegen references.
orbital_reference_test = executable(
    'kumo-orbital-reference-test',
    'tests/orbital_reference_test.cpp',
    dependencies: kumo_dep,
)
foreach orbital : [
    'o2_alpha_mo_3',
    'o2_alpha_mo_9',
    'o2_alpha_mo_10',
    'o2_beta_mo_3',
    'o2_beta_mo_7',
    'o2_beta_mo_8',
    'acrolein_mo_13',
    'acrolein_mo_15',
    'acrolein_mo_16',
    'acrolein_mo_18',
]
    parts = orbital.split('_')
    molecule = parts[0]
    channel = molecule == 'o2' ? parts[1] : 'alpha'
    index = parts[parts.length() - 1]
    test(
        'orbital-reference-' + orbital,
        orbital_reference_test,
        args: [meson.current_source_dir() / (
            'tests/data/' + molecule + '.fchk.gz'
        ), meson.current_source_dir() / (
            'tests/data/mo_reference_cubes/' + orbital + '.cube'
        ), channel, index],
        timeout: 300,
    )
endforeach

kumo_cli_test = find_program('tests/cli_list_test.sh')
test(
    'cli-list-restricted',
    kumo_cli_test,
    args: [kumo_cli, meson.current_source_dir() / 'tests/data/acrolein.fchk.gz', 'restricted'],
)
test(
    'cli-list-unrestricted',
    kumo_cli_test,
    args: [kumo_cli, meson.current_source_dir() / 'tests/data/o2.fchk.gz', 'alpha'],
)

# Exercise the installed command-line interface on a real FCHK file.
kumo_cli_e2e = find_program('tests/cli_end_to_end.sh')
test(
    'cli-end-to-end',
    kumo_cli_e2e,
    args: [kumo_cli, meson.current_source_dir() / 'tests/data/o2.fchk.gz', meson.project_version()],
    timeout: 300,
)

# Progress callbacks: monotonic, documented totals, identical results.
test(
    'progress',
    executable('kumo-progress-test', 'tests/progress_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Derived density properties: ELF, reduced gradient, sign(lambda2) rho.
test(
    'derived-properties',
    executable(
        'kumo-derived-properties-test',
        'tests/derived_properties_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/hybrids/h2o_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Derived density properties against Multiwfn reference Cube files.
test(
    'multiwfn-reference',
    executable(
        'kumo-multiwfn-reference-test',
        'tests/multiwfn_reference_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/hybrids/h2o_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/multiwfn',
    ],
)

# The AO shell cutoff zeroes only AOs and derivatives below ao_cutoff.
ao_shell_cutoff_test_exe = executable(
    'kumo-ao-shell-cutoff-test',
    'tests/ao_shell_cutoff_test.cpp',
    dependencies: kumo_dep,
)
test(
    'ao-shell-cutoff',
    ao_shell_cutoff_test_exe,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 120,
)
# The large HBC case is in the 'slow' suite, which the sanitizer CI job skips.
test(
    'ao-shell-cutoff-hbc',
    ao_shell_cutoff_test_exe,
    args: [meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz'],
    suite: 'slow',
    timeout: 120,
)

# Batched AO derivatives agree with the single-point evaluator.
test(
    'ao-derivatives-batch',
    executable(
        'kumo-ao-derivatives-batch-test',
        'tests/ao_derivatives_batch_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
        meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz',
    ],
)

# Spatial AO derivatives (analytic implementation versus numerical oracle).
test(
    'ao-spatial-derivatives',
    executable(
        'kumo-ao-spatial-derivatives-test',
        'tests/ao_derivatives_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Analytic spatial derivatives of signed molecular orbitals.
test(
    'mo-spatial-derivatives',
    executable(
        'kumo-mo-spatial-derivatives-test',
        'tests/mo_derivatives_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Analytic spatial derivatives of spin-resolved electron densities.
test(
    'density-spatial-derivatives',
    executable(
        'kumo-density-spatial-derivatives-test',
        'tests/density_derivatives_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Independent MO-occupation reconstruction and blocked-versus-pointwise derivatives.
test(
    'density-mo-block-consistency',
    executable(
        'kumo-density-mo-block-consistency-test',
        'tests/density_mo_consistency_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Optional profiling of spatial derivatives; not part of the CI test suite.
if benchmarks_enabled
    executable(
        'kumo-spatial-benchmark',
        'benchmarks/spatial_derivatives_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

# Gzip regression tests use committed compressed references and materialize
# plain-text copies only inside a temporary test directory.
test(
    'fchk-gzip',
    executable('kumo-fchk-gzip-test', 'tests/fchk_gzip_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
        meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_631gdp.fchk.gz',
        meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_def2svp.fchk.gz',
        meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_def2tzvp.fchk.gz',
    ],
)

# Compare density grids with 1/2/4 OpenMP workers (also runs in serial builds).
parallel_density_test = executable(
    'kumo-parallel-density-test',
    'tests/parallel_density_test.cpp',
    dependencies: kumo_dep,
)
test(
    'parallel-density',
    parallel_density_test,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)
parallel_derivatives_test = executable(
    'kumo-parallel-derivatives-test',
    'tests/parallel_derivatives_test.cpp',
    dependencies: kumo_dep,
)
test(
    'parallel-derivatives',
    parallel_derivatives_test,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Evaluator point path: bit-identical results across OpenMP worker counts.
evaluator_parallel_test = executable(
    'kumo-evaluator-parallel-test',
    'tests/evaluator_parallel_test.cpp',
    dependencies: kumo_dep,
)
test(
    'evaluator-parallel',
    evaluator_parallel_test,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# MO grid and spatial-derivative consistency across OpenMP worker counts.
parallel_orbitals_test = executable(
    'kumo-parallel-orbitals-test',
    'tests/parallel_orbitals_test.cpp',
    dependencies: kumo_dep,
)
test(
    'parallel-orbitals',
    parallel_orbitals_test,
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# AO-P versus AO-MO derivative benchmark.
if benchmarks_enabled
    executable(
        'kumo-derivative-backend-comparison',
        'benchmarks/derivative_backend_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

test(
    'evaluator-api',
    executable('kumo-evaluator-api-test', 'tests/evaluator_api_test.cpp', dependencies: kumo_dep),
    args: [meson.current_source_dir() / 'tests/data/o2.fchk.gz'],
)

test(
    'grid-validation',
    executable(
        'kumo-grid-validation-test',
        'tests/grid_validation_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [meson.current_source_dir() / 'tests/data/o2.fchk.gz'],
)

# Fused/individual evaluator requests against validated 1.3 kernels.
test(
    'evaluator-fusion',
    executable(
        'kumo-evaluator-fusion-test',
        'tests/evaluator_fusion_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

test(
    'parallel-blocks',
    executable(
        'kumo-parallel-blocks-test',
        ['tests/parallel_blocks_test.cpp'],
        include_directories: include_directories('src'),
        dependencies: kumo_dep,
    ),
)

# Boys function against quadrature and closed forms.
test(
    'boys-function',
    executable(
        'kumo-boys-function-test',
        ['tests/boys_function_test.cpp'],
        include_directories: include_directories('src'),
        dependencies: kumo_dep,
    ),
)

# Electrostatic potential against closed forms, the multipole expansion and
# the Poisson equation.
electrostatic_potential_exe = executable(
    'kumo-electrostatic-potential-test',
    ['tests/electrostatic_potential_test.cpp'],
    include_directories: include_directories('src'),
    dependencies: kumo_dep,
)
test(
    'electrostatic-potential',
    electrostatic_potential_exe,
    args: [
        meson.current_source_dir() / 'tests/data/hybrids/h2o_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 120,
)
# HBC 6-31G(d,p) adds Cartesian d and SP shells. It is in the 'slow' suite,
# which the sanitizer CI job skips: it took 90 s under AddressSanitizer
# locally, several times more on the hosted runner.
test(
    'electrostatic-potential-hbc',
    electrostatic_potential_exe,
    args: [meson.current_source_dir() / 'tests/data/hbc/hbc_b3lyp_631gdp.fchk.gz'],
    suite: 'slow',
    timeout: 120,
)

# Internal dependency planning, without changing numerical kernels.
test(
    'evaluation-plan',
    executable(
        'kumo-evaluation-plan-test',
        ['tests/evaluation_plan_test.cpp'],
        include_directories: include_directories('src'),
        dependencies: kumo_dep,
    ),
)

# Direct numerical regression against the legacy AO-P derivative kernel.
test(
    'shared-ao-validation',
    executable(
        'kumo-shared-ao-validation-test',
        'tests/shared_ao_validation_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Optional single-thread algorithmic comparison; run manually, not in CI.
if benchmarks_enabled
    executable(
        'kumo-shared-ao-benchmark',
        'benchmarks/shared_ao_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

# First-order AO path against full second-order contractions.
test(
    'selective-ao-derivatives',
    executable(
        'kumo-selective-ao-derivatives-test',
        'tests/selective_ao_derivatives_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Benchmark selective first-order versus second-order AO derivatives.
if benchmarks_enabled
    executable(
        'kumo-selective-ao-benchmark',
        'benchmarks/selective_ao_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

# Validate AO-MO eligibility without changing production contractions.
test(
    'backend-selection',
    executable(
        'kumo-backend-selection-test',
        'tests/backend_selection_test.cpp',
        include_directories: include_directories('src'),
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
)

# Experimental AO-MO kernel, independent of Evaluator selection.
test(
    'ao-mo-kernel',
    executable('kumo-ao-mo-kernel-test', 'tests/ao_mo_kernel_test.cpp', dependencies: kumo_dep),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Manual three-way AO-MO numerical diagnostic.
if benchmarks_enabled
    executable('kumo-ao-mo-diagnostic', 'tests/ao_mo_diagnostic.cpp', dependencies: kumo_dep)
endif

# Explicit AO_P / AO_MO Evaluator backend selection.
test(
    'evaluator-backend',
    executable(
        'kumo-evaluator-backend-test',
        'tests/evaluator_backend_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Manual end-to-end AO_P vs AO_MO comparison through Evaluator.
if benchmarks_enabled
    executable(
        'kumo-evaluator-backend-benchmark',
        'benchmarks/evaluator_backend_benchmark.cpp',
        dependencies: kumo_dep,
    )
endif

# Deterministic, conservative automatic backend selection.
test(
    'automatic-backend',
    executable(
        'kumo-automatic-backend-test',
        'tests/automatic_backend_test.cpp',
        include_directories: include_directories('src'),
        dependencies: kumo_dep,
    ),
)

# Measure AUTO decisions against both explicit backends.
if benchmarks_enabled
    executable(
        'kumo-automatic-backend-benchmark',
        'benchmarks/automatic_backend_benchmark.cpp',
        include_directories: include_directories('src'),
        dependencies: kumo_dep,
    )
endif

# Analytic AO overlaps and Gaussian MO metric regression.
test(
    'analytic-overlap',
    executable(
        'kumo-analytic-overlap-test',
        'tests/analytic_overlap_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Independent Gauss-Hermite primitive-product overlap reference.
test(
    'analytic-overlap-reference',
    executable(
        'kumo-analytic-overlap-reference-test',
        'tests/analytic_overlap_reference_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Expose the stored spin-resolved FCHK 1-RDM without duplicating state.
test(
    'one-particle-density',
    executable(
        'kumo-one-particle-density-test',
        'tests/one_particle_density_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# S-metric subspace projectors (factorized representation).
test(
    'subspace-projector',
    executable(
        'kumo-subspace-projector-test',
        'tests/subspace_projector_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Metric projector stress tests.
test(
    'subspace-projector-stress',
    executable(
        'kumo-subspace-projector-stress-test',
        'tests/subspace_projector_stress_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Quadrature-weighted AO overlap accumulation (no grid policy).
test(
    'weighted-overlap',
    executable(
        'kumo-weighted-overlap-test',
        'tests/weighted_overlap_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Diagnostic convergence of spatial quadrature toward analytic AO overlap.
test(
    'weighted-overlap-quadrature',
    executable(
        'kumo-weighted-overlap-quadrature-test',
        'tests/weighted_overlap_quadrature_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Independent atomic quadrature validation.
test(
    'atomic-radial-angular-quadrature',
    executable(
        'kumo-atomic-radial-angular-quadrature-test',
        'tests/atomic_radial_angular_quadrature_test.cpp',
        dependencies: kumo_dep,
    ),
    timeout: 300,
)

# Higher-order single-centre quadrature convergence.
test(
    'atomic-radial-angular-convergence',
    executable(
        'kumo-atomic-radial-angular-convergence-test',
        'tests/atomic_radial_angular_convergence_test.cpp',
        dependencies: kumo_dep,
    ),
    timeout: 300,
)

# End-to-end single-center AO overlap quadrature.
test(
    'weighted-overlap-single-center',
    executable(
        'kumo-weighted-overlap-single-center-test',
        'tests/weighted_overlap_single_center_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [
        meson.current_source_dir() / 'tests/data/o2.fchk.gz',
        meson.current_source_dir() / 'tests/data/acrolein.fchk.gz',
    ],
    timeout: 300,
)

# Becke multicenter overlap diagnostic.
test(
    'weighted-overlap-multicenter',
    executable(
        'kumo-weighted-overlap-multicenter-test',
        'tests/weighted_overlap_multicenter_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [meson.current_source_dir() / 'tests/data/o2.fchk.gz'],
    timeout: 300,
)

# Generic affine-grid trilinear field sampler.
test(
    'field-sampler',
    executable('kumo-field-sampler-test', 'tests/field_sampler_test.cpp', dependencies: kumo_dep),
)

# Validate grid geometry independently of physical scale.
test(
    'grid-geometry',
    executable('kumo-grid-geometry-test', 'tests/grid_geometry_test.cpp', dependencies: kumo_dep),
    timeout: 120,
)

# Mutable field samples preserve geometry and sampler consistency.
test(
    'field-mutability',
    executable(
        'kumo-field-mutability-test',
        'tests/field_mutability_test.cpp',
        dependencies: kumo_dep,
    ),
)

# Preserve distinct atomic number and explicit center charge in Cube output.
test(
    'atom-charge-cube',
    executable(
        'kumo-atom-charge-cube-test',
        'tests/atom_charge_cube_test.cpp',
        dependencies: kumo_dep,
    ),
    timeout: 120,
)

# General scalar Cube I/O and metadata round-trip.
test(
    'cube-io',
    executable('kumo-cube-io-test', 'tests/cube_io_test.cpp', dependencies: kumo_dep),
    timeout: 120,
)

# Read independently supplied Gaussian Cube reference files.
test(
    'cube-real-reference',
    executable(
        'kumo-cube-real-reference-test',
        'tests/cube_real_reference_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [meson.current_source_dir() / 'tests/data'],
    timeout: 120,
)

# Generic scalar-field integration, independent of trilinear interpolation.
test(
    'field-integrals',
    executable('kumo-field-integrals-test', 'tests/field_integrals_test.cpp', dependencies: kumo_dep),
)

# Analytic AO Cartesian moments: s/p/d, pure shells, translations, and FCHK ordering.
test(
    'cartesian-ao-moments',
    executable(
        'kumo-cartesian-ao-moments-test',
        'tests/cartesian_ao_moments_test.cpp',
        dependencies: kumo_dep,
    ),
    args: [meson.current_source_dir() / 'tests/data/o2.fchk.gz'],
    timeout: 120,
)

# Independent molecular dipole benchmark from Gaussian FCHK records.
# Reference values come from Gaussian, not from a second KUMO integration.
gaussian_dipole_reference = executable(
    'kumo-gaussian-dipole-reference-test',
    'tests/gaussian_dipole_reference_test.cpp',
    dependencies: kumo_dep,
)
test(
    'gaussian-dipole-reference',
    gaussian_dipole_reference,
    args: [
        meson.current_source_dir() / 'tests/data/hybrids/h2o_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/hybrids/nh3_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/hybrids/ch4_b3lyp_def2tzvp.fchk.gz',
    ],
    timeout: 120,
)

# Independent Gaussian FCHK quadrupole reference; test all six tensor components.
gaussian_quadrupole_reference = executable(
    'kumo-gaussian-quadrupole-reference-test',
    'tests/gaussian_quadrupole_reference_test.cpp',
    dependencies: kumo_dep,
)
test(
    'gaussian-quadrupole-reference',
    gaussian_quadrupole_reference,
    args: [
        meson.current_source_dir() / 'tests/data/hybrids/h2o_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/hybrids/nh3_b3lyp_def2tzvp.fchk.gz',
        meson.current_source_dir() / 'tests/data/hybrids/ch4_b3lyp_def2tzvp.fchk.gz',
    ],
    timeout: 120,
)
