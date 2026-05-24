# CLAUDE.md — GRAB project conventions

## Release status — pre-1.0, no backward compatibility constraint

GRAB has not yet been publicly released; there are no external users and no
deployed binaries in the field.  Therefore **do not** preserve backward
compatibility when redesigning interfaces:

- CLI flags, value enumerations, default behaviors, and file formats may be
  renamed, restructured, or removed outright when a cleaner design is
  available.
- Do not add aliases, deprecation warnings, fallback parsers, or "legacy
  mode" branches for flags or values that have been renamed.  Let the old
  spelling fall through to the dispatcher's generic "unknown option" error.
  Do not add bespoke `--old-flag → --new-flag` redirect branches either:
  pre-1.0 there are no external scripts to redirect, and the redirect
  branches just turn into clutter the next time the flag is renamed.
- Do not preserve old function signatures, struct layouts, or output column
  orderings "in case downstream consumers depend on them" — there are none.
- The only durability constraint is reproducibility within a single
  release: a given GRAB binary, given the same inputs and flags, must
  produce the same outputs across runs.

Apply this principle when weighing design changes: choose the cleaner long-
term interface and break the older one, rather than carrying parallel paths.

## Build & dependency model

GRAB is a **pure C++17 application**. The repository is fully self-contained:

- All third-party libraries are **vendored** under `third_party/` (Eigen, Boost
  headers-only subset, zlib, zstd, libdeflate, plink2/pgenlib, bgen, htslib).
- The build pulls **nothing** from the system: no `apt-get`, no `brew install`,
  no `vcpkg`, no Conan. The only requirements are a C++17 compiler (g++ /
  clang++ / MinGW g++) and `make`.
- `Makefile` produces a single statically-linked binary (`build/grab2`) that runs
  on **Linux, macOS, and Windows (MSYS2/MinGW)** with no shared-library
  dependencies. Distribution = copy the binary.

**Therefore:** any change you make must keep this property. When considering
performance, parallelism, or new features:

- **Do not** add OpenMP (`-fopenmp`, `#pragma omp parallel`, `#pragma omp simd`).
  OpenMP requires `libomp`/`libgomp`, which is a runtime dependency that
  Apple clang on macOS does not ship by default.
- **Do not** add Intel TBB, oneAPI, MKL, OpenBLAS, or other parallel/numerical
  runtimes.
- **Do not** introduce `find_package`, CMake-based system probes, or any
  mechanism that consults `/usr/lib`, `/usr/local`, or Homebrew paths.
- **Do not** add new third-party libraries unless you also vendor them under
  `third_party/` and confirm they build on all three platforms with the
  existing toolchain.

**Do** use:

- C++17 standard library (`<thread>`, `<atomic>`, `<mutex>`,
  `<condition_variable>`, `<future>`) for parallelism.
- The vendored Eigen for vectorized math — Eigen ships its own SIMD intrinsics
  (SSE/AVX on x86, NEON on ARM) and does not require any external runtime.
- The existing **runtime SIMD dispatch** infrastructure in
  `src/util/simd_dispatch.hpp` and `src/util/simd_math.hpp`. New SIMD kernels
  should follow the SPAsqr pattern (`src/spasqr/spasqr.cpp`):
  - Write `_avx512` and `_avx2` variants tagged with
    `__attribute__((target("avx2,fma")))` /
    `__attribute__((target("avx2,avx512f,avx512vl,fma")))`.
  - Keep a scalar fallback that compiles on any architecture (including ARM).
  - Resolve via `simdLevel()` at first call.
  - x86 detection guarded by `#if defined(__x86_64__) || defined(_M_X64)`.

**AVX2 fallback is mandatory for every AVX-512 kernel.** GRAB targets
x86-64-v3 hardware as the supported deployment baseline; AVX-512 is a
performance enhancement on top of that, not a deployment requirement.
Therefore:

- Every code path that contains an `_avx512` variant must also ship an
  `_avx2` variant.  The runtime dispatcher in `simd_dispatch.hpp` picks
  `_avx2` whenever `__builtin_cpu_supports("avx512...")` returns false,
  which is the common case on consumer Intel chips and on many cloud
  instances.
- Do not write AVX-512-only kernels and rely on the scalar fallback as
  the secondary path.  On hosts without AVX-512 that would silently lose
  the SIMD speed-up that the AVX2 variant provides.
- When you add a new SIMD kernel, the review checklist is: scalar +
  `_avx2` + `_avx512` variants, plus a `simdLevel()` dispatch site.
  Anything missing the AVX2 tier is a bug, not an optimization
  opportunity for later.

## Architecture matrix the build supports

| Platform | x86_64                      | arm64 (Apple Silicon, ARM Linux) |
| -------- | --------------------------- | -------------------------------- |
| Linux    | AVX-512 / AVX2 / scalar     | NEON via Eigen + scalar          |
| macOS    | AVX-512 / AVX2 / scalar     | NEON via Eigen + scalar          |
| Windows  | AVX-512 / AVX2 / scalar     | (untested)                       |

The Makefile auto-detects platform (`uname -s`), architecture (`uname -m`),
and AVX2 availability. `GRAB_MARCH=-march=native` is the default; override
with `GRAB_MARCH=-march=x86-64-v2` for portable distribution binaries.

## Source layout

- `src/cli/`     — argument parsing, dispatch, help.
- `src/engine/`  — marker / LOCO chunk-level work-stealing thread pool
                   (uses `std::thread`, no external runtime).
- `src/io/`      — sparse GRM, subject data, file readers.
- `src/geno_factory/` — genotype format readers (plink, pgen, bgen, vcf).
- `src/{spacox, spagrm, spamix, spasqr, wtcoxg, localplus}/` — methods.
- `src/util/`    — math helpers, SIMD dispatch, vectorized exp/log,
                   logging, text scanners, IQR outlier detection
                   (`outlier.hpp` shared by spamix / wtcoxg).

## Conventions for new code

- Eigen `VectorXd` / `ArrayXd` for vector math; prefer Eigen array ops
  (`(a.array() * b.array()).sum()`) over hand-written loops — Eigen's
  expression templates dispatch to SIMD automatically.
- For per-variant hot loops with transcendental functions (exp/log), follow
  the SPAsqr SIMD-dispatch pattern; do not assume the compiler will
  auto-vectorize through scalar `std::exp` calls.
- Chunk-level parallelism is provided by the marker engine
  (`src/engine/marker.cpp`) via `std::thread` and an atomic chunk counter —
  per-method code does not spawn threads.

## Building, testing, packaging

- `make -j$(nproc)` — builds `build/grab2`.
- `make clean` — removes `build/`.
- The binary is the deliverable. There is no install step, no shared library,
  no headers exposed to users. Users download or build the binary and run it.

## Example scripts — `examples/tutorial.sh` and `examples/baseline.sh`

The repository ships two end-to-end example scripts with distinct
audiences and obligations.

### `examples/tutorial.sh` — user-facing walkthrough

`tutorial.sh` is a minimal, copy-pasteable demonstration of each
analysis method.  Every command lists only the mandatory inputs plus a
phenotype list and an output prefix; all other knobs fall back to
`grab2`'s built-in defaults.  Two phenotypes are exercised per method
to keep the output tables small while still covering both quantitative
and survival code paths.  This script is **not** a regression baseline;
its output is not pinned to a hash.

### `examples/baseline.sh` — exhaustive regression script

`baseline.sh` exercises every utility mode (`--cal-af-coef`,
`--cal-pairwise-ibd`, `--int-pheno`) and every analysis method
(SPACox, SPAmix, SPAGRM, SAGELD, SPAsqr in both **score** and **wald**
modes, WtCoxG, LEAF) against the bundled `examples/1kg.*` fixtures,
writing all artifacts under `examples_output/`.

The script has two purposes that must both be preserved when it is
edited:

1. **Documentary.**  Every command spells out every command-line flag
   that the method accepts, with each numeric or categorical knob set to
   its built-in default value.  Reading the script tells a new
   contributor exactly which flags are available and what each one
   defaults to when omitted.  When a new flag is added to a method or
   utility, the corresponding block in `examples/baseline.sh` must gain
   a line listing the flag at its default value; when a flag is removed
   or renamed, the line must be updated accordingly.

2. **Regression baseline.**  After any refactor — shared engine, SIMD
   kernels, null-model fitting, genotype readers, output formatting, or
   per-method code — re-run `examples/baseline.sh` and confirm that the
   resulting `examples_output/*` artifacts are byte-identical (or
   numerically identical up to documented tolerance) to the
   pre-refactor baseline.  A passing build is not sufficient evidence
   that a refactor preserved behavior; output equivalence is.

The compression codec varies across blocks by design, so that a single
pass through the script exercises all three output-writer paths:
plain text (SPACox, SPAsqr-wald), gzip (WtCoxG, LEAF), and zstd
(everything else).  Do not collapse the codec to a single setting when
editing the script.

The `--int-pheno` block sits between SAGELD and SPAsqr: it produces
`examples_output/1kg.int.txt`, an inverse-normal-transformed phenotype
file containing only the `Quantitative` and `Time` columns.  SPAsqr
consumes this file via `--pheno` and pulls the remaining covariates
(`MALE`, `PC1..PC4`) from the original phenotype file via `--covar`;
this exercises the disjoint-pheno/covar loading path in `SubjectData`.

The SPAsqr **wald** block restricts to 100 variants via
`--extract examples/spasqr_wald_extract`, because the per-marker QR
refit is appreciably slower than score mode.  The 100-line ID file is
checked into the repository under `examples/` to keep the regression
result reproducible.

The cross-format SPAGRM block at the bottom of the script converts the
bundled `.pgen` fixture to BED, BCF, and BGEN with `plink2 --make-bed
/ --export bcf / --export bgen-1.2` and runs SPAGRM on each input with
the same `--extract` / `--exclude` / `--keep` / `--remove` filter
lists, then asserts byte-identity across the four readers in
`src/geno_factory/`.  This serves as both the cross-reader regression
and the regression for the shared `geno_factory::filterMarkersByIds`
ID-filter helper.

## Shared engine code is validated — do not modify when debugging other methods

SPAsqr, SPAmix, SPACox, SPAGRM, SAGELD, WtCoxG, and LEAF are currently
passing end-to-end tests.  Together they exercise every facet of the shared
engine infrastructure: the fused-GEMM path (SPAsqr, SPAGRM, SPAmix,
SAGELD-pheno mode, WtCoxG, LEAF), the `MissBatch` non-fuseable path
(SPACox), the single-phenotype engine (SAGELD residual mode), the LOCO
engine (SPAsqr-LOCO), the multi-phenotype engine (SPAsqr, SPAGRM, SPAmix,
SPACox), the per-cluster sub-method pattern (LEAF), and the runtime
SIMD-dispatch pattern.  Their collective success is a strong signal that
**the common code below is correct**.

The only method currently under active debugging is **SPAmixLocalPlus**.
When debugging it, do not suspect or modify the shared infrastructure
listed below — the bug is almost certainly in the method-specific code,
not the shared engine:

- `src/engine/marker.cpp`, `src/engine/marker.hpp`, `src/engine/marker_impl.hpp`
  — `markerEngine`, `multiPhenoEngine`, `multiPhenoEngineRange`,
  chunk-level work-stealing thread pool, fused union-level GEMM
  (`AugResid^T × GBatch_union`), `FusedStatsGroup` QC sharing,
  `MissBatch` extraction path for non-fuseable phenotypes.
- `src/engine/loco.cpp`, `src/engine/loco.hpp` — `locoEngine`,
  Regenie / LDAK-KVIK `.loco` parsers, per-chromosome task rebuild loop.
- `src/util/simd_dispatch.hpp`, `src/util/simd_math.hpp` — runtime
  AVX2 / AVX-512 dispatch and vectorized exp/log kernels.
- `src/util/null_model.{hpp,cpp}` — `parseRegressionModel`, the unified
  null-model fitting engine driving the `--pheno-name + --regression-model`
  path for the seven validated methods.
- `src/geno_factory/` — genotype decoding (plink / pgen / bgen / vcf):
  SPAsqr and SPAmix exercise all four readers through the same
  `GenoCursor` interface.
- The `MethodBase` interface contract in `src/engine/marker.hpp`
  (`clone`, `prepareChunk`, `getResultVec`, `getResultBatch`,
  `supportsFusedGemm`, `fillUnionResiduals`, `fillResidualSums`,
  `processScoreBatch`).

If a method under debug misbehaves, look first at its own per-method file
(score centering, null-model fitting, residual construction, p-value
computation, output formatting, QC thresholds it sets itself).  Changing
the shared engine to "fix" a method bug will break the seven validated
methods and is the wrong direction; if a shared-engine change is genuinely
required, re-run the regression tests for SPAsqr, SPAmix, SPACox, SPAGRM,
SAGELD, WtCoxG, and LEAF before committing.
