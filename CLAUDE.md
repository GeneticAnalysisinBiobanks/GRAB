# CLAUDE.md — GRAB project conventions

## Release status — pre-1.0, no backward compatibility constraint

GRAB has not yet been publicly released; there are no external users and no
deployed binaries in the field.  Therefore **do not** preserve backward
compatibility when redesigning interfaces:

- CLI flags, value enumerations, default behaviors, and file formats may be
  renamed, restructured, or removed outright when a cleaner design is
  available.
- Do not add aliases, deprecation warnings, fallback parsers, or "legacy
  mode" branches for flags or values that have been renamed.  Reject the old
  spelling with an unambiguous error message that points to the new one.
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
- `Makefile` produces a single statically-linked binary (`build/grab`) that runs
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

- `make -j$(nproc)` — builds `build/grab`.
- `make clean` — removes `build/`.
- The binary is the deliverable. There is no install step, no shared library,
  no headers exposed to users. Users download or build the binary and run it.

## Shared engine code is validated — do not modify when debugging other methods

SPAsqr is currently passing end-to-end tests. Because SPAsqr is the most
demanding consumer of the shared engine infrastructure (it exercises the
fused-GEMM path, the multi-phenotype engine, LOCO, and the SIMD-dispatch
pattern), a working SPAsqr is a strong signal that **the common code below
is correct**. When debugging any other method (SPAGRM / SPACox /
WtCoxG / SPAmix / LEAF), do not suspect or modify these files — the bug is
almost certainly in the method-specific code, not the shared infrastructure:

- `src/engine/marker.cpp`, `src/engine/marker.hpp`, `src/engine/marker_impl.hpp`
  — `markerEngine`, `multiPhenoEngine`, `multiPhenoEngineRange`,
  chunk-level work-stealing thread pool, fused union-level GEMM
  (`AugResid^T × GBatch_union`), `FusedStatsGroup` QC sharing,
  `MissBatch` extraction path for non-fuseable phenotypes.
- `src/engine/loco.cpp`, `src/engine/loco.hpp` — `locoEngine`,
  Regenie / LDAK-KVIK `.loco` parsers, per-chromosome task rebuild loop.
- `src/util/simd_dispatch.hpp`, `src/util/simd_math.hpp` — runtime
  AVX2 / AVX-512 dispatch and vectorized exp/log kernels.
- `src/geno_factory/` — genotype decoding (plink / pgen / bgen / vcf):
  SPAsqr drives all four readers through the same `GenoCursor` interface.
- The `MethodBase` interface contract in `src/engine/marker.hpp`
  (`clone`, `prepareChunk`, `getResultVec`, `getResultBatch`,
  `supportsFusedGemm`, `fillUnionResiduals`, `fillResidualSums`,
  `processScoreBatch`).

If a non-SPAsqr method misbehaves, look first at its own per-method file
(score centering, null-model fitting, residual construction, p-value
computation, output formatting, QC thresholds it sets itself). Changing
the shared engine to "fix" a method bug will break SPAsqr and is the wrong
direction; if a shared-engine change is genuinely required, re-run SPAsqr
regression tests before committing.
