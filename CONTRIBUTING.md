# Contributing

## Maintain the rules

- Maintain development policy in this file; link to it from other instruction files.
- Keep build recipes in `Makefile`, release jobs in `.github/workflows/release.yml`,
  command-line definitions in `src/cli/flags.hpp`, and the version in `VERSION`.
- Update each rule at its authoritative location; do not copy it into another guide.
- Write concise imperatives. Omit rationale, implementation history, and progress logs.

## Build

- Use C++17 and the existing Makefile on Linux, macOS, and Windows with MSYS2/MinGW.
- Vendor application dependencies under `third_party/`; do not add system package
  discovery, network downloads, OpenMP, or external numerical runtimes to the build.
- Validate new vendored dependencies on all supported platforms.
- Preserve compiler and flag overrides: `CXX`, `CC`, `CPPFLAGS`, `CXXFLAGS`,
  `CFLAGS`, and `LDFLAGS`.
- Use the same Makefile for local and release builds. Follow the release workflow's
  platform settings when validating release builds.
- Preserve native CPU tuning for local x86 builds, the AVX2 release baseline,
  and ARM support. Keep release runtime linkage settings platform-specific.
- Run commands from the repository root. Set a job count appropriate for the host:

  ```bash
  make -j8
  ```

- Run `make clean` before changing compilers or build flags; preserve any needed
  artifacts outside `build/` and `dist/` first.

## Implement changes

- Replace obsolete alpha interfaces directly; do not add compatibility aliases,
  deprecated paths, fallback parsers, or duplicate signatures.
- Implement each shared numerical procedure once in `src/util/`; pass method-specific
  inputs through adapters.
- Fix method-local defects in method code. Treat required shared-code changes as
  explicit work affecting every consumer; validate every affected method.
- Include SPAmix when validating SPAmixPlus changes and all callers when changing
  shared adapters.
- Use Eigen arrays and vectors for vector arithmetic.
- Use C++17 threading facilities. Keep per-marker parallelism in the marker engine;
  do not spawn threads inside method marker loops.
- Make results reproducible across repeated runs, thread counts, and chunk sizes.
- For single instruction, multiple data (SIMD) kernels, provide scalar, AVX2, and
  AVX-512 implementations and select them through `simdLevel()`.
- Guard x86 intrinsics by architecture and use the existing function target
  attributes in `src/util/spa_cgf.cpp`; retain portable scalar compilation.
- Use `src/util/simd_math.hpp` for vectorized exponential and logarithm operations.
- Keep numerical tolerances consistent with their quantities' units. Preserve
  p-values under phenotype rescaling.
- Compute output p-values in the log domain as `LOG10P` (negative base-10 logarithm
  of the p-value); use the shared distribution, combination, and inversion helpers.
- Use the status definitions in `src/util/spa.hpp` and the output contract in
  `src/cli/flags.hpp`. Report unavailable results as `NA`; do not silently clamp,
  saturate, or replace failed results.
- Use only documented fallback estimators with an explicit status and validated
  accuracy. Do not substitute an estimator when no valid statistic exists.
- Keep SPAmixLocalPlus formatting consistent with `src/engine/marker_impl.hpp`.

## Select checks

Configure output paths once per shell session:

```bash
export TEMP="$PWD/tests_output"
mkdir -p "$TEMP"
export TMPDIR="$TEMP"
GRAB_CHECK_DIR=$(mktemp -d "$TEMP/grab2-checks.XXXXXX")
set -o pipefail
```

- Keep executables and build dependencies in `build/`.
- Write test data, logs, reports, and other disposable artifacts under `tests_output/`.
- Keep `examples/baseline.sh` and `examples/tutorial.sh` outputs under
  `examples_output/baseline/` and `examples_output/tutorial/`, respectively.
- Run standalone test executables from the repository root; use
  `tests/test_output.hpp` for their default fixture directories under `tests_output/`.
- Use the run directory printed by `make test` or `make bench` to inspect fixtures
  and per-executable logs; retain it until review is complete.

Run every applicable row before submitting changes.

| Change | Required checks |
| --- | --- |
| Documentation only | Check paths, links, commands, and agreement with the implementation. |
| C++ source or headers | Build and run `make test`. |
| Analysis methods, shared engine, numerical code, readers, or output | Also run the baseline comparison below. |
| Output-column selection | Also run `tests/output_columns.py`. |
| Performance-sensitive code | Also compare relevant benchmarks before and after. |
| Statistical calculations | Also check independent reference values, scale invariance, and null calibration. |
| Build or dependency configuration | Also validate local and release builds on affected platforms. |

- Add a focused regression test for a bug fix where practical; confirm that it
  fails before the fix and passes afterward.
- Name C++ tests `tests/*_test.cpp` and benchmarks `tests/bench_*.cpp`.
- Declare additional test object dependencies as `TESTOBJS_<name>` in `Makefile`.
- Keep assertions enabled in test builds.
- Run the C++ suite separately from the example scripts:

  ```bash
  make -j8 test 2>&1 | tee "$GRAB_CHECK_DIR/test.log"
  ```

## Compare complete workflows

- Provide Bash, Python 3, `plink2`, `zstd`/`zstdcat`, and the shell utilities used
  by `examples/baseline.sh`, including `md5sum`.
- Start each baseline run with a fresh `examples_output/baseline/` directory;
  archive existing results before running the script.
- Before editing code, build and save reference results:

  ```bash
  make -j8
  bash examples/baseline.sh
  mv examples_output/baseline examples_output/baseline.before
  ```

- After editing code, rebuild, test, and compare:

  ```bash
  make -j8
  make -j8 test 2>&1 | tee "$GRAB_CHECK_DIR/test-after.log"
  bash examples/baseline.sh
  bash tests/regress.sh examples_output/baseline.before examples_output/baseline \
    2>&1 | tee "$GRAB_CHECK_DIR/regress.log"
  ```

- Use an unused reference directory; retain its revision and build settings.
- If no reference exists, generate it from the previous revision in a separate
  checkout using the same toolchain and settings.
- Inspect every baseline `[FAIL]` message; do not accept the script's exit status
  alone. Review changes in columns, rows, missing values, signs, and numeric values.
- Preserve results by default. Require byte identity or agreement within a
  documented tolerance for behavior-preserving changes; use `--rtol 1e-9` unless
  the change specifies a different validated tolerance.
- For intentional numerical changes, document the expected differences, check
  independent reference values, and validate every affected method before accepting
  a new baseline. Do not add compatibility paths or relax tolerances to hide failures.
- Review reported changes in missing values and `LOG10P` magnitudes explicitly;
  do not treat a successful comparison exit status as their acceptance criterion.
- Reject negative or infinite `LOG10P` values and inconsistent result/status pairs.
- Account separately for `converted/1kg.log` timestamps and resource estimates.
  Leave the comparison tool's binary variant file date normalization unchanged.
- Compare `fit.ibd.zst` normally between revisions; exclude it only from checks
  across thread counts until its existing reduction-order variation is fixed.

## Run additional checks

- Run output-selection checks in a new scratch directory:

  ```bash
  python3 tests/output_columns.py --work-dir "$GRAB_CHECK_DIR/output-columns" \
    2>&1 | tee "$GRAB_CHECK_DIR/output-columns.log"
  ```

- Run benchmarks before and after the change on the same machine with the same
  build settings and comparable load:

  ```bash
  make -j8 bench 2>&1 | tee "$GRAB_CHECK_DIR/bench.log"
  ```

- Use `tests/make_synthetic.py` for genotype cohorts, `tests/make_lanc_null.py`
  for local-ancestry null cohorts, and `tests/make_fallback.py` with
  `tests/fallback_report.py` for fallback checks. Follow each script's usage and
  provide NumPy for the generators.
- Set generator destinations and save console reports as follows:

  | Script | Output argument | Console report |
  | --- | --- | --- |
  | `tests/make_synthetic.py` | `--out "$GRAB_CHECK_DIR/synthetic"` | `$GRAB_CHECK_DIR/synthetic.log` |
  | `tests/make_fallback.py` | `--out "$GRAB_CHECK_DIR/fallback"` | `$GRAB_CHECK_DIR/fallback.log` |
  | `tests/make_lanc_null.py` | `"$GRAB_CHECK_DIR/lanc-null"` | `$GRAB_CHECK_DIR/lanc-null.log` |
  | `tests/fallback_report.py` | Pass result tables under `tests_output/`. | `$GRAB_CHECK_DIR/fallback-report.log` |

- Append `2>&1 | tee REPORT` to each command, replacing `REPORT` with its console
  report path. Put analysis outputs from generated cohorts under `tests_output/` with `--out`.
- Add a suitable longitudinal null cohort when changing SAGELD's statistical
  calculations; do not claim null calibration from the bundled examples alone.

## Maintain examples

- Keep `examples/tutorial.sh` minimal: mandatory inputs, phenotype selections,
  and output prefixes.
- List accepted options at their defaults in `examples/baseline.sh`; omit boolean
  switches whose default is off. Update affected blocks when flags change.
- Preserve plain-text, gzip, and zstd coverage, separate phenotype/covariate inputs,
  the checked-in SPAsqr Wald subset, and comparisons across all four genotype readers.

## Submit changes

- Record the behavior change, affected methods, commands run, results, and any
  checks not completed in the change description.
