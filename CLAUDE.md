# CLAUDE.md — GRAB project conventions

## Release status — alpha, no backward compatibility constraint

GRAB is approaching its 2.0-alpha release.  There are no production
deployments and no external scripts depending on the current interface, so
**do not** add back-compat scaffolding when redesigning user-visible
surfaces:

- CLI flags, value enumerations, default behaviors, and file formats may be
  renamed, restructured, or removed outright when a cleaner design emerges.
- Do not add aliases, deprecation warnings, fallback parsers, or per-flag
  redirect branches.  Let the old spelling fall through to the dispatcher's
  generic "unknown option" error — redirects just turn into clutter the
  next time the flag is renamed.
- Do not preserve old function signatures, struct layouts, or output column
  orderings "in case downstream consumers depend on them" — there are none.
- The only durability constraint is reproducibility within a single
  release: a given GRAB binary, given the same inputs and flags, must
  produce the same outputs across runs.

Apply this principle when weighing design changes: choose the cleaner long-
term interface and break the older one, rather than carrying parallel paths.

## The unification principle — one implementation per shared procedure

`feat/spa-unify` is a **refactor branch**, and the refactor is not confined
to the saddlepoint.  Its governing rule is:

> When two or more methods perform the same statistical or numerical
> procedure, the repository must contain exactly **one** implementation of
> it, placed in `src/util/`, with each method reduced to supplying its own
> inputs and consuming the shared result.

The saddlepoint tier (`src/util/spa.hpp`, `src/util/spa_cgf.{hpp,cpp}`) is
the completed instance of the rule: six hand-copied root finders, five
hand-copied binomial CGFs and several hand-copied tail assemblies collapsed
into one of each, with the per-method code reduced to a tier-3 CGF adapter.
The same rule governs every other procedure shared across methods — p-value
representation and the columns that carry it, statistic-from-p-value
inversion, p-value combination (Cauchy, meta, conditional), null-model
fitting, Wald refits, QC thresholds, output-column contracts — whether or
not that procedure has been unified yet.  A duplicate discovered during
this work is a defect to be removed, not a parallel path to be maintained.

The **currently active** instance of the rule is the p-value representation
itself: `dev-notes/methods/log10p_unify/`.  It makes `LOG10P` (= −log10 P)
the sole p-value column, deletes the linear `P` column and the parallel
linear tail assembly that produces it, moves every p-value inversion,
combination and pooling into the log domain, and re-partitions
`SPA_STATUS`.  While that project is in flight, the "Saddlepoint output
columns" section below describes the *outgoing* contract; the incoming one
is `dev-notes/methods/log10p_unify/00_decisions.md`, and where the two
disagree the plan wins.

Unification is also the performance strategy, not merely a tidiness
argument.  One implementation is the only one that can afford a
cancellation-free formulation, a scalar/AVX2/AVX-512 triple, and a test
suite proportionate to its blast radius; N copies get none of those,
because no single copy justifies the investment.

### Changed numbers are the expected outcome, not a regression

Replacing N per-method copies with one shared implementation **changes
results in every method that had a copy**, and it must: the copies did not
agree with each other, so no single implementation can reproduce all of
them.  Each copy also carried its own accumulated defects — a different
convergence tolerance, a cancellation-prone cumulant algebra, a guard at a
different threshold, a `std::min(1.0, NaN)` that laundered a saddlepoint
failure into a perfectly null marker.  Correcting those moves p-values by
design.

Therefore, on this branch:

- **Do not** read a diff in `examples_output/` as evidence of a defect.  A
  unification stage that produced no diff at all would be the surprising
  outcome and would itself warrant investigation.
- **Do not** add a compatibility switch, a "legacy" code path, a
  per-method constant, or a tolerance fudge in order to make a previous
  number reappear.  That recreates exactly the duplication the branch
  exists to remove, and it is forbidden by the alpha-release rule above.
- **Do** be able to name, for every column that moved, the unification
  decision that moved it and the direction the change should have gone.
  An **unexplained** diff is the failure signal; a diff as such is not.

### What must be preserved instead of the numbers

Because the numbers are free to move, the invariants carry the weight.  A
unification stage is acceptable when all of the following still hold, and
each of them is checkable:

1. **Statistical validity.**  Under the null the method's p-values remain
   uniform, and its type-I error at genome-wide thresholds is unchanged.
   `tests/make_synthetic.py` and `tests/make_lanc_null.py` generate the
   pure-null-by-construction cohorts this is measured on; the bundled
   `examples/1kg.*` fixture is far too small to resolve calibration.
2. **Direction of correctness.**  Where old and new differ materially, the
   new value is the one that agrees with an independent reference — a
   higher-precision evaluation, an analytic special case, or a
   `tests/spa_reference.hpp` implementation written for clarity rather
   than speed.  "It changed" is not a finding; "it changed and the new
   value is the correct one, by this reference" is.
3. **Dimensional soundness and scale equivariance.**  Rescaling a
   phenotype is a change of units and must not move a p-value at all; see
   `src/util/spa.hpp` below.  `scale_equivariance_*` in
   `tests/spa_solver_test.cpp` pins this as bit-identity.
4. **Failure is named, never fabricated.**  A procedure that cannot
   produce an answer reports `NA` together with a status naming the
   reason.  It never substitutes a plausible number, and it never floors,
   clamps or saturates a quantity into a value a reader would mistake for
   a measurement.  This is the single most damaging class of defect the
   repository has had to remove, and every unification stage must leave it
   removed.

   What this rule forbids is the **silent** substitution.  Substituting a
   different, well-defined, documented estimator and naming the
   substitution in a status column is a different act and is permitted:
   the `log10p_unify` project reports the normal-approximation tail when
   the saddlepoint fails, under `SPA_STATUS` codes 3–6 which exist for no
   other purpose, and the reader recovers the unsubstituted subset with
   `SPA_STATUS <= 2`.  The boundary is that a fallback needs (a) a genuine
   estimator behind it, (b) a status value that names it, and (c) a
   documented warning where the substituted estimator is known to be
   biased — for the normal fallback, that it is anti-conservative at low
   MAC, which is exactly where the saddlepoint fails.  A fallback lacking
   any of the three is fabrication.  Note also that "no statistic exists"
   (no informative subject, monomorphic stratum, `Var(S) <= 0`) is never
   eligible for a fallback: there is nothing to fall back *to*, and those
   rows stay `NA`.
5. **Reproducibility.**  A given binary, given the same inputs and flags,
   produces the same outputs across runs, thread counts and chunk sizes.

`tests/regress.sh` is the instrument for reading a diff, and on this branch
it is diagnostic rather than pass/fail: read its per-column `max |abs|`,
`max rel` and `max dlog10P` report and ask whether every moving column is
accounted for by (1)-(5), not whether the report came back empty.

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

**AVX2 fallback is mandatory for every AVX-512 kernel.** GRAB's supported
deployment baseline is x86-64-v3; AVX-512 is a performance enhancement on
top of that, not a deployment requirement.  The runtime dispatcher in
`simd_dispatch.hpp` picks `_avx2` whenever AVX-512 is unavailable (the
common case on consumer Intel chips and many cloud instances), so an
AVX-512-only kernel with a scalar fallback silently loses the SIMD speed-
up on those hosts.  Every new SIMD kernel must ship the full scalar +
`_avx2` + `_avx512` triple plus a `simdLevel()` dispatch site; anything
missing the AVX2 tier is a bug, not a follow-up.

## Two distribution paths share one Makefile

A single `Makefile` serves both source-compile users and the GitHub
Actions release-build workflow.  The two paths differ only in a small
number of variable values; do not introduce a second Makefile.

- **Source-compile users (`make -j`).**  Default `GRAB_MARCH := -march=native`
  tunes the resulting binary to the build host's CPU, enabling Eigen's
  and the compiler's AVX-512 auto-vectorization on capable hardware in
  addition to the runtime-dispatched SPAsqr/SPAmix/SPAGRM kernels.  This
  is the optimal path for users who run GRAB on the same machine they
  build it on.
- **GitHub Actions release builds
  (`make -j GRAB_MARCH=-march=haswell`).**  Pins the baseline ISA at
  AVX2/FMA/BMI2 so the published binary runs on any x86-64 machine with
  a 2013-or-newer CPU.  `haswell` is used in place of the PSABI-level
  name `x86-64-v3` because the manylinux2014 container ships GCC 10,
  which predates `x86-64-v3` (introduced in GCC 11); the two flags
  target the same ISA so codegen is equivalent.  The runtime SIMD
  dispatcher in `simd_dispatch.hpp` still picks the AVX-512 variants
  of hand-written kernels on capable hosts (because
  `__attribute__((target(...)))` emits all variants regardless of
  `-march`); only auto-vectorized code is capped at AVX2.

Both `CXX`/`CC` and the standard `CPPFLAGS`/`CXXFLAGS`/`CFLAGS`/`LDFLAGS`
are honored via `?=` and trailing-append patterns, so an external build
environment can inject sysroot-aware compilers and hardening flags
(`-fstack-protector-strong`, `-D_FORTIFY_SOURCE=2`,
`-fdebug-prefix-map=…`, rpath additions, etc.) without further
intervention.  A representative release-build invocation, as used by
[.github/workflows/release.yml](.github/workflows/release.yml):

```bash
make -j$(nproc) \
    GRAB_MARCH="-march=haswell" \
    SIMD_FLAGS="-mavx2 -mbmi -mbmi2 -mlzcnt -mfma" \
    STATIC_LIBS="-static-libstdc++ -static-libgcc"
make install PREFIX="${PREFIX}"
```

The Linux release build runs inside the `quay.io/pypa/manylinux2014_*`
container so that the resulting binary's glibc baseline is 2.17 (RHEL 7
era), maximizing HPC-cluster compatibility.  `STATIC_LIBS` is set
exactly because the manylinux container ships a newer libstdc++ than
the target hosts; static-linking libstdc++/libgcc cuts the runtime
dependency back to glibc only.

For `linux-aarch64` and `osx-arm64`, omit `SIMD_FLAGS` and let
`GRAB_MARCH` default to empty — Eigen's NEON path is selected
automatically and the scalar fallback covers everything else.  For
`macos-*`, omit `STATIC_LIBS` because Apple clang does not accept
`-static-libstdc++` (libc++ is provided by the system).

When renaming or removing a Makefile variable, audit both paths: keep
the source-compile default produce-the-fastest behavior intact, and
keep the release workflow minimal (no per-platform branches inside the
workflow that the Makefile could absorb).

## Architecture matrix the build supports

| Platform | x86_64                      | arm64 (Apple Silicon, ARM Linux) |
| -------- | --------------------------- | -------------------------------- |
| Linux    | AVX-512 / AVX2 / scalar     | NEON via Eigen + scalar          |
| macOS    | AVX-512 / AVX2 / scalar     | NEON via Eigen + scalar          |
| Windows  | AVX-512 / AVX2 / scalar     | (untested)                       |

The Makefile auto-detects platform (`uname -s`), architecture (`uname -m`),
and AVX2 availability. `GRAB_MARCH=-march=native` is the default; override
with `GRAB_MARCH=-march=x86-64-v3` for portable distribution binaries, or
`-march=haswell` on GCC 10 and older, which is the spelling the release
workflow uses.  Do not drop to `x86-64-v2`: it predates AVX2 and is below
the deployment baseline stated above, so the auto-vectorized and
Eigen-expressed code would lose AVX2 while the runtime-dispatched kernels
kept it.

## Source layout

- `src/cli/`     — argument parsing, dispatch, help.
- `src/engine/`  — marker / LOCO chunk-level work-stealing thread pool
                   (uses `std::thread`, no external runtime).
- `src/io/`      — sparse GRM, subject data, file readers.
- `src/geno_factory/` — genotype format readers (plink, pgen, bgen, vcf).
- `src/{spacox, spagrm, spagxe, spamix, spasqr, wtcoxg, localplus}/` —
                   methods.  Two methods have no directory of their own:
                   **LEAF** is `src/wtcoxg/leaf.{cpp,hpp}` (a per-cluster
                   composition of WtCoxG) and **SAGELD** is
                   `src/spagrm/sageld*.{cpp,hpp}` (it reuses the SPAGRM
                   null model).  `src/localplus/` additionally holds the
                   `.lanc` local-ancestry format (`lanc_io`,
                   `lanc_convert_rfmix`) behind `--make-lanc` / `--lanc`.
- `src/util/`    — the shared numerical and infrastructural tier: math
                   helpers and distribution wrappers (`math_helper`),
                   the saddlepoint tier (`spa.hpp`, `spa_cgf`),
                   fixed-effects meta-analysis (`meta_pvalue.hpp`, used
                   by LEAF), standard-model Wald tests (`wald`, used by
                   the SPAGxE Branch-B leg), unified null-model fitting
                   (`null_model`), least squares (`regression`), the INT
                   transform (`int_pheno`), read-only mmap (`mmap_ro`),
                   compressed text I/O (`text_stream`), SIMD dispatch and
                   vectorized exp/log, logging, text scanners, and IQR
                   outlier detection (`outlier.hpp`, shared by spamix,
                   spagxe, wtcoxg/leaf and localplus).
  - `src/util/spa.hpp` — the shared saddlepoint tier: the bracketed and
    safeguarded root solver, the Barndorff-Nielsen tail kernel (linear and
    log domain), the two-sided assembly, and the `Status`/`SPA_STATUS`
    encoding (see "Saddlepoint output columns" below).

    **Every threshold in the solver is dimensionally sound, and this is a
    pinned property, not a stylistic preference.** Multiplying a phenotype —
    or a `--resid-name` residual vector — by a constant `c` is a pure change
    of units: the statistic becomes `c*S`, the saddlepoint root becomes
    `zeta/c`, and `w`, `v`, `r*` and the p-value are all exactly invariant.
    Residual tolerances are therefore built only from quantities carrying the
    units of `s` (`sqrt(K'')`, `\|s\|`, `K''*\|t\|`) and abscissa distances only
    from quantities carrying the units of `t` (`\|t\|`, `\|K'/K''\|`). An
    absolute constant in either place makes the reported p-value depend on
    the units the phenotype happens to be recorded in; the predecessor's
    `stepTol * max(1, \|t\|)` bracket tolerance and `bracketStep * max(1,
    \|init\|)` coarse probe were exactly that, and moved 44 % of SPA-branch
    p-values under a `x1e8` rescaling of a heavy-tailed phenotype.
    `scale_equivariance_*` in `tests/spa_solver_test.cpp` pins the property
    as **bit-identity** of `P`, `LOG10P` and `SPA_STATUS` under a
    power-of-two rescaling, which is exact in binary floating point. Do not
    introduce a new absolute constant into the solver without checking it
    against that test.
  - `src/util/spa_cgf.{hpp,cpp}` — the three shared binomial CGF variants
    (`binomUniform`, `binomIndiv`, `binomHapcount`), each with the mandated
    scalar + AVX2 + AVX-512 triple and a `simdLevel()` dispatch site.
  - Five tier-3 adapters sit on top of `spa_cgf`:
    `src/{spacox,spagrm,spamix,wtcoxg}/*_cgf.hpp` and
    `src/localplus/spamixlocalp_cgf.hpp`. SPAsqr consumes
    `spa_cgf::binomUniform` directly, with no adapter file of its own.
    The adapters are **not** one-per-method: `spamix_cgf.hpp` serves
    SPAmix, SPAmixPlus, SPAGxE and SPAGxEmix, and `wtcoxg_cgf.hpp` serves
    WtCoxG and LEAF.  Editing either touches more than one method, which
    matters because the debugging advice below directs you to "the
    method's own adapter" — for four of the ten entry points that file is
    shared.

## Saddlepoint output columns — `LOG10P` and `SPA_STATUS`

> **Under active revision.**  `dev-notes/methods/log10p_unify/` deletes the
> `P` column described below, makes `LOG10P` the sole p-value column across
> all ten method entry points, and replaces the seven-value `spa::Status`
> with a nine-value partition that separates "the saddlepoint failed and
> the normal tail was reported instead" (codes 3–6) from "no test exists"
> (code 8) and "a step downstream of the saddlepoint failed" (code 7).
> Read `00_decisions.md` and `02_inventory.md` there before touching any
> output column.  Everything below this box describes the contract as it
> stands *before* that project lands.

Every method that *reports a saddlepoint p-value* emits a matching
`LOG10P` and `SPA_STATUS` beside it.  The plain spelling `P LOG10P Z
Z_Norm BETA SE SPA_STATUS` is used by SPACox, SPAGRM and
SPAmix/SPAmixPlus; the rest prefix or suffix the two names to match the
p-value they qualify — `*_EXT` / `*_NOEXT` in WtCoxG, `meta_*` and
`cl<N>_*` in LEAF, `*_tau<t>` in SPAsqr, `anc<K>_*` in the unpinned
SPAmixLocalPlus, `*_Gx<E>` in SPAGxE/SPAGxEmix.  The two columns are:

- `LOG10P` — `-log10(P)`, carried through the tail assembly in the log
  domain (`spa::bnTailLog` + `spa::combineTailsLog`, via `math::pnormLog`)
  so it stays meaningful past the point where the linear-scale tail
  underflows to exactly zero.
- `SPA_STATUS` — `static_cast<uint8_t>(spa::Status)` (`src/util/spa.hpp`),
  naming the outcome of the saddlepoint that produced the row:

  | Value | Token | Meaning |
  | ----- | ----- | ------- |
  | 0 | `OK` | Converged: the residual `\|K'(zeta) - s\|` met the solver's stated criterion at the returned root, re-verified there against the terminal cumulants; no guard fired. |
  | 1 | `MAXITER` | The solver stopped without meeting that criterion — the iteration or bracket-expansion budget ran out, or the bracket closed to the last representable digits of its own endpoints with the residual still above tolerance. Both mean no root was located to the requested accuracy. |
  | 2 | `GUARD_TEMP` | `zeta*s - K(zeta) < 0`, so `w` is not real. |
  | 3 | `GUARD_CURV` | `K''(zeta) <= 0`, so `v` is not real. |
  | 4 | `GUARD_W` | `\|w\|` at or below the removable-singularity threshold (`spa::kWSingularity`, `1e-3`) — a **degraded success**: `Phi(+/-w)` is returned in place of the `r*` correction, and `P` is a number. |
  | 5 | `NONFINITE` | `zeta`, a cumulant, or `r*` left the reals. |
  | 6 | `NORMAL` | `\|Z\| <= --spa-z-threshold`; the saddlepoint was never attempted and `P` is the plain normal-approximation tail. |

  `P` and `LOG10P` are `NA` for every status other than 0, 4 and 6 — a
  saddlepoint failure is reported as a named reason rather than silently
  substituted with a fabricated p-value. `SPA_STATUS` is emitted as this
  integer code rather than the string token `spa::statusName` produces,
  because `MethodBase::getResultVec` / `getResultBatch` hand the marker
  engine a `std::vector<double>` per row with no string channel
  (`src/engine/marker_impl.hpp`'s `numToChars` formats every result cell),
  so a non-numeric result column is structurally impossible rather than
  merely unconventional.  Introducing one would be an `src/engine/`
  redesign — permissible as a unification change with the full
  re-validation that implies, never as a side effect of method work.

### The contract is not yet universal — known gaps

The rule above describes the columns the saddlepoint tier itself produces.
Several p-values in the tree still escape it, and closing the gap is
unification work, not an optional cleanup.  **All six gaps below are
scheduled to close in `dev-notes/methods/log10p_unify/`** — see
`04_validation.md` §6 for the stage that closes each, and for the two that
project deliberately leaves open.  Do not open a separate effort against
any of them.

- **SAGELD** emits `P_G Z_G BETA_G SE_G` and the matching `_Gx<E>` group
  (`src/spagrm/sageld.cpp`) with **no** `LOG10P` and **no** `SPA_STATUS`,
  even though its p-value comes from `SPAGRMClass::getMarkerPvalFromScore`,
  which returns all three.  The two columns are discarded at the call site.
- **`P_CCT`** (SPAsqr score and wald, SPAGxE, SPAGxEmix) has no `LOG10P`
  companion, because `math::cauchyCombine` is evaluated entirely on the
  linear scale.
- **`cl<i>_P_BAT`** in LEAF (`src/wtcoxg/leaf.cpp`) has no `LOG10P`
  companion.
- **`LOG10P_EXT` / `LOG10P_NOEXT`** in WtCoxG are computed as
  `-std::log10(p)` from the already-assembled linear ratio in
  `wtcoxg_cond::conditionalP` (`src/wtcoxg/conditional_p.hpp`), not in the
  log domain, so they become `+Inf` rather than a magnitude once the
  conditional p underflows.
- The **`Z` column** is derived from the linear `P` through
  `math::zFromPval`, whose `math::qnorm` clamps its argument at `1e-300`.
  It therefore saturates at `|Z| = 37.047096` for every marker at or below
  that p, while the adjacent `LOG10P` keeps rising — the two columns
  disagree in exactly the regime where the result matters most.
  `math::qchisq` carries the same clamp and saturates at `1373.87`, which
  bounds the per-cluster weight in `math::metaPvalueScorePool` and the
  variance recovered from a p-value in WtCoxG's conditional branches.

Until these are unified, do not describe the `P` / `LOG10P` / `Z` /
`SPA_STATUS` contract to a user as if it held everywhere.

## Conventions for new code

- Eigen `VectorXd` / `ArrayXd` for vector math; prefer Eigen array ops
  (`(a.array() * b.array()).sum()`) over hand-written loops — Eigen's
  expression templates dispatch to SIMD automatically.
- For per-variant hot loops with transcendental functions (exp/log), follow
  the SPAsqr SIMD-dispatch pattern; do not assume the compiler will
  auto-vectorize through scalar `std::exp` calls.
- Per-marker and per-chunk parallelism is owned by the marker engine
  (`src/engine/marker.cpp`) via `std::thread` and an atomic chunk counter.
  Per-method code must not spawn threads **inside the per-marker loop**.
  Null-model fitting and preprocessing phases do run their own pools
  (`util/null_model.cpp`, `spagrm/{grm_null,ibd}.cpp`,
  `spamix/indiv_af.cpp`, `spacox/spacox_cgf.hpp`, `spasqr/spasqr.cpp`,
  `wtcoxg/{wtcoxg,leaf}.cpp`, `spagrm/sageld.cpp`), and SPAmixLocalPlus
  bypasses the marker engine entirely with a pool of its own.  When adding
  one, respect the same reproducibility constraint the engine does: the
  result must not depend on thread count or scheduling order.

## Building, testing, packaging

- `make -j$(nproc)` — builds `build/grab2`.
- `make clean` — removes `build/`.
- The binary is the deliverable. There is no install step, no shared library,
  no headers exposed to users. Users download or build the binary and run it.
- `make test` and `make bench` are the developer targets. `make test` builds
  and runs every `tests/*_test.cpp` binary, running **all** of them and
  exiting non-zero at the end if any failed (it does not stop at the first
  failure); `make bench` builds and runs every `tests/bench_*.cpp` binary.
  Both are `.PHONY`: `tests/` sits outside the source tree that feeds
  `build/grab2`, so nothing under it is linked into the shipped binary and
  the "binary is the deliverable" property above is unaffected. There is no
  third-party test framework — `tests/tinytest.hpp` is a single
  self-contained assertion header — and tests build with the same
  `GRAB_CXXFLAGS` (`-march`, SIMD flags, optimization level) as `src/`, with
  `-DNDEBUG` filtered back out so assertions fire. As of the saddlepoint
  unification, `make test` runs ten binaries: seven SPA suites
  (`spa_solver_test`, `spa_cgf_test`, `spacox_cgf_test`, `spagrm_cgf_test`,
  `spamix_cgf_test`, `wtcoxg_cgf_test`, `spamixlocalp_cgf_test`) plus the
  three `.lanc` suites (`lanc_simd_test`, `lanc_roundtrip_test`,
  `lanc_convert_rfmix_smoke_test`); `make bench` currently builds one
  (`bench_spa_cgf`).  Two coverage gaps are known rather than accidental:
  `src/wtcoxg/conditional_p.hpp` is exercised only indirectly through
  `wtcoxg_cgf_test`, and there is no SPAGxE suite at all even though
  SPAGxE/SPAGxEmix are pinned in `examples/baseline.sh`.
- `tests/make_synthetic.py` and `tests/make_lanc_null.py` generate the
  synthetic null cohorts the calibration gates use — pure-null-by-construction
  genotype/phenotype/GRM fixtures at a scale (50 000 subjects x 20 000
  markers at the script's defaults, larger via `--n` / `--m`) the bundled
  `examples/1kg.*` fixture is too small to resolve kernel-level or
  calibration-level change on. They are invoked directly when reproducing a
  calibration; they are not part of `make test`.

## Example scripts — `examples/tutorial.sh` and `examples/baseline.sh`

The repository ships two end-to-end example scripts with distinct
audiences and obligations.

### `examples/tutorial.sh` — user-facing walkthrough

A minimal, copy-pasteable demonstration of each analysis method.  Every
command lists only the mandatory inputs, a phenotype list, and an output
prefix; all other knobs fall back to defaults.  Two phenotypes per method
keep tables small while covering both quantitative and survival code
paths.  Output is **not** pinned to a hash.

### `examples/baseline.sh` — exhaustive regression script

Exercises every utility mode (`--cal-af-coef`, `--cal-pairwise-ibd`,
`--int-pheno`) and every validated analysis method against the bundled
`examples/1kg.*` fixtures, writing all artifacts under
`examples_output/baseline/`: SPACox in fit, residual and longitudinal
modes; SPAmix with pre-computed and on-the-fly AF; SPAGRM; SAGELD in fit,
residual and GALLOP modes; SPAsqr in both **score** and **wald** modes;
WtCoxG; LEAF; SPAGxE base and sparse-GRM; SPAGxEmix — plus per-chromosome
LOCO for SPACox, SPAmix, SPAGRM, SPAsqr and WtCoxG.  Two purposes must be
preserved when editing it:

1. **Documentary.**  Every command spells out every flag the method
   accepts, with each numeric or categorical knob set to its built-in
   default.  When a flag is added, the corresponding block must gain a
   line listing it at its default; when a flag is renamed or removed,
   the line must be updated.

2. **Change ledger.**  After any change — shared engine, shared numerical
   tier, SIMD kernels, null-model fitting, genotype readers, output
   formatting, or per-method code — re-run the script and diff the result
   against the previous tree.  A passing build is not evidence about
   behavior.  Note the obligation is to **account for** the diff, not to
   have none: per "Changed numbers are the expected outcome" above, a
   unification stage is supposed to move results, and the deliverable is
   a per-column explanation, not a zero.

`tests/regress.sh examples_output.base examples_output` (backed by
`tests/regress.py`) is the instrument: for every artifact that differs it
reports added/removed/reordered columns, the per-column maximum absolute
and relative difference, the maximum change in -log10(P) ("max dlog10P"),
and NA/exactly-zero/sign-change counts by direction, rather than merely
flagging that two files disagree.  Plain `cmp` is not the right tool for
this repository and never was: a stricter convergence tolerance, a
cancellation-free cumulant algebra, a repaired guard or a newly shared
procedure each move p-values below `cmp`'s bit-for-bit bar, so `cmp`
reports "different" for every file and tells you nothing about which
change did it.

Know what the script actually gates on, because on this branch it will
usually exit non-zero and that is not by itself a finding:

- The **only** numeric gate is `max rel > --rtol`, with `--rtol`
  defaulting to `1e-9` and applied uniformly to every numeric column — `P`,
  `LOG10P`, `Z`, `BETA` and `SE` alike.
- `max dlog10P` is computed and printed but **never** compared to a
  threshold; there is no flag to set one.  The per-stage budget it is read
  against lives in `dev-notes/methods/spa_unify/03_stages.md`, not in the
  script.
- **Structural findings always fail**: a file present on one side only, a
  differing row count, an added/removed/reordered column, or any change in
  a non-numeric column.  Every stage of this branch that introduces a
  `LOG10P` or `SPA_STATUS` column therefore exits non-zero by construction.
- Changes in missingness (`NA` gained or lost) are counted and printed but
  do **not** fail on their own, even though a marker crossing into or out
  of `NA` is a behavioral change and must be named in the ledger.
- `baseline/converted/1kg.log` must be discounted **by hand**: it is
  plink2's own log, embedding a wall-clock timestamp and a memory
  estimate, and is not reproducible run to run even on an unmodified tree.
  `tests/regress.py` has no exclusion list, so delete or ignore that file
  before reading the report.

Then read the report column by column.  The questions it has to answer
are: which columns moved, by how much on the `-log10(P)` scale a GWAS
result is actually read on, whether any cell changed to or from `NA`,
whether any cell changed sign, and whether any column appeared or
disappeared (a contract change, which must be intentional and documented
in `src/cli/flags.hpp`).  A `max dlog10P` of 1e-9 across a whole method is
rounding.  A handful of markers moving by 0.5 with the rest at 1e-12 is a
guard or branch threshold that changed, and the ledger must say which one
and why the new side is right.

The compression codec varies across blocks by design, so a single pass
exercises all three output-writer paths: plain text (the SPACox fit,
residual and longitudinal blocks), gzip (SAGELD-GALLOP, WtCoxG,
WtCoxG-LOCO, LEAF, SPAGxE sparse-GRM), and zstd (everything else,
including SPAsqr in **both** score and wald modes).  Do not collapse the
codec to a single setting.

Three blocks deserve specific attention when editing:

- **`--int-pheno` → SPAsqr.**  The `--int-pheno` block produces
  `examples_output/baseline/int_pheno.txt` (an inverse-normal-transformed
  file with only `Quantitative` and `Time` columns).  SPAsqr then consumes
  it via `--pheno` while pulling the remaining covariates (`MALE`,
  `PC1..PC4`) from the original phenotype file via `--covar`, which
  exercises the disjoint-pheno/covar loading path in `SubjectData`.
- **SPAsqr wald.**  Restricted to 8 variants via
  `--extract examples/spasqr_wald_extract`, since the per-marker QR
  refit is appreciably slower than score mode.  The 8-line ID file
  is checked in to keep the regression reproducible.
- **Cross-format SPAGRM.**  The trailing block converts the bundled
  `.pgen` fixture to BED, BCF, and BGEN (via `plink2 --make-bed /
  --export bcf / --export bgen-1.2`), runs SPAGRM on each input with
  identical `--extract` / `--exclude` / `--keep` / `--remove` filters,
  and asserts byte-identity across the four readers in
  `src/geno_factory/`.  This regresses both the readers and the shared
  `geno_factory::filterMarkersByIds` helper.

## Shared code is changed deliberately, never incidentally

This section is not in tension with the unification principle above; the
two are the same rule seen from opposite ends.  Shared code **is** the
object of this branch's work, and rewriting it is the point.  What is
forbidden is changing it *incidentally* — as a side effect of chasing one
method's bug — because every such change silently alters the nine other
methods that depend on it, and none of them are in front of you at that
moment.

The distinction is the entry condition:

- **A unification change** starts from the shared tier, states which
  methods it will move and in which direction, changes the shared code
  once, and re-validates **every** consumer against the invariants in
  "What must be preserved" before it is committed.  This is the normal
  mode of work on `feat/spa-unify`.
- **A debugging change** starts from one method behaving wrongly.  Here
  the shared tier is almost certainly not the culprit, and editing it is
  the wrong direction.

SPAsqr, SPAmix, SPACox, SPAGRM, SAGELD, WtCoxG, LEAF, SPAGxE and
SPAGxEmix all currently pass end-to-end and are pinned in
`examples/baseline.sh`.  Together they exercise every facet of the shared
engine infrastructure:

| Facet | Methods that exercise it |
| ----- | ------------------------ |
| fused union-level GEMM | SPAsqr, SPAGRM, WtCoxG, LEAF, and SPAmix in pre-computed-AF mode only |
| non-fused `getResultBatch` / `MissBatch` | SPACox, SAGELD, SPAGxE/SPAGxEmix, and SPAmix in on-the-fly-AF mode |
| LOCO engine | SPACox, SPAGRM, SPAmix, SPAsqr (both modes), WtCoxG |
| multi-phenotype engine | SPACox, SPAGRM, SPAmix, SPAsqr, SAGELD-pheno mode |
| single-phenotype engine | SAGELD residual mode |
| per-cluster sub-method pattern | LEAF |
| per-marker Wald refit | SPAsqr-wald, SAGELD, SPAGxE |
| runtime SIMD dispatch | SPAsqr, SPAmix, SPAGRM, and the shared `spa_cgf` kernels |

Nine methods agreeing is strong evidence that the common code below is
correct, and correspondingly strong evidence that a change made to it
while debugging one of them will break the other eight.

**Ten method entry points** depend on the saddlepoint tier
(`src/util/spa.hpp`, `src/util/spa_cgf.{hpp,cpp}`, added by the `spa:`
commits on `feat/spa-unify`): SPACox, SPAGRM, SAGELD, SPAsqr, SPAmix,
SPAmixPlus, SPAGxE, SPAGxEmix, WtCoxG, LEAF, and the unpinned
SPAmixLocalPlus — a larger and different set than the engine roster above,
because SAGELD reaches the tier through `SPAGRMClass`, SPAGxE/SPAGxEmix
through `spamix_cgf.hpp`, and LEAF through `wtcoxg_cgf.hpp`.  A
method-specific saddlepoint bug is almost certainly in that method's
tier-3 CGF adapter — the code supplying `K`, `K'` and `K''` — not in the
shared solver, tail kernel, or CGF variants.  Note that **SAGELD discards
`LOG10P` and `SPA_STATUS`** at the call site (`src/spagrm/sageld.cpp`), so
a tier change can move SAGELD p-values with no status column present to
diagnose it with.

The methods currently under active development are **SPAmixPlus** and
**SPAmixLocalPlus**.  Both remain callable via `--method spamixplus` /
`--method spamixlocalplus` but are hidden from `grab2 --help`; neither
appears in `examples/baseline.sh`, and their outputs are not pinned.

SPAmixPlus is **not** a separate implementation: `--method spamix` and
`--method spamixplus` share `src/spamix/spamixplus.{cpp,hpp}` and the same
`SPAmixPlusMethod` class, differing only at the dispatch branch.  Any edit
made while debugging SPAmixPlus changes the **pinned SPAmix output** as
well and must be re-gated with `tests/regress.sh`.  SPAmixLocalPlus is
genuinely separate — it does not go through the marker engine at all — but
it does keep a hand-mirrored second copy of the result-formatting policy
(`src/localplus/spamixlocalp.cpp`), which must be kept in step with
`src/engine/marker_impl.hpp`.  When debugging either method, do not
suspect or modify the shared infrastructure listed below — the bug is
almost certainly in the method-specific code:

- `src/engine/marker.cpp`, `src/engine/marker.hpp`, `src/engine/marker_impl.hpp`
  — `markerEngine`, `multiPhenoEngine`, `multiPhenoEngineRange`,
  chunk-level work-stealing thread pool, fused union-level GEMM
  (`AugResid^T × GBatch_union`), `FusedStatsGroup` QC sharing,
  `MissBatch` extraction path for non-fuseable phenotypes.
- `src/engine/loco.cpp`, `src/engine/loco.hpp` — `locoEngine`,
  Regenie / LDAK-KVIK `.loco` parsers, per-chromosome task rebuild loop.
- `src/util/simd_dispatch.hpp`, `src/util/simd_math.hpp` — runtime
  AVX2 / AVX-512 dispatch and vectorized exp/log kernels.
- `src/util/spa.hpp` — the shared saddlepoint tier: `solveSaddlepoint` (the
  bracketed-and-safeguarded Newton solver), `bnTailLog` (the
  Barndorff-Nielsen modified-signed-root tail, log domain), `combineTailsLog`
  and `assemble` (the two-sided assembly and the decision-D5 fallback), and
  the `Status`/`SPA_STATUS` encoding described above.  The linear twins
  `bnTail`, `combineTails` and `normalTwoSided` were deleted by
  `log10p_unify` Stage 3; a method that still emits a `P` column derives it
  as `spa::pFromNegLog10P` of the assembled `-log10(P)`, and that helper goes
  with the column in Stage 8.
- `src/util/spa_cgf.{hpp,cpp}` — the three shared binomial CGF variants
  (`binomUniform`, `binomIndiv`, `binomHapcount`), each with the mandated
  scalar + AVX2 + AVX-512 triple and a `simdLevel()` dispatch site.
- `src/util/math_helper.{hpp,cpp}` — the distribution tier every method
  inverts or evaluates through: `pnorm`, `pnormLog`, `qnorm`, `qchisq`,
  `pt`, `zFromPval`, `cauchyCombine`, `pmvnorm2dHalfRect`.
- `src/util/meta_pvalue.hpp` — `metaPvalueScorePool`, the inverse-variance
  score pooling LEAF combines its per-cluster results with.
- `src/util/wald.{hpp,cpp}` — the shared Wald refit fitters used by
  SPAsqr-wald, SAGELD and SPAGxE.
- `src/util/null_model.{hpp,cpp}` — `parseRegressionModel`, the unified
  null-model fitting engine driving the `--pheno-name + --regression-model`
  path for SPACox, SPAGRM, SPAmix/SPAmixPlus, SPAGxE/SPAGxEmix, WtCoxG,
  LEAF and SPAmixLocalPlus.  SPAsqr fits its own smoothed-QR null model
  (`src/spasqr/qmme.{cpp,hpp}`) and SAGELD its own longitudinal one
  (`src/spagrm/sageld_fit.{cpp,hpp}`); neither accepts
  `--regression-model`.
- `src/geno_factory/` — genotype decoding (plink / pgen / bgen / vcf):
  SPAsqr and SPAmix exercise all four readers through the same
  `GenoCursor` interface.
- The `MethodBase` interface contract in `src/engine/marker.hpp`
  (`clone`, `resultSize`, `getHeaderColumns`, `prepareChunk`,
  `getResultVec`, `getResultBatch`, `preferredBatchSize`,
  `supportsFusedGemm`, `fusedGemmColumns`, `fillUnionResiduals`,
  `fillResidualSums`, `processScoreBatch`).  `getHeaderColumns` is the
  hook the output-column contract above rests on.

If a method under debug misbehaves, look first at its own per-method file
(score centering, residual construction, the tier-3 CGF adapter, p-value
assembly, output formatting, QC thresholds it sets itself).  Reaching into
the shared tier to "fix" one method's bug will break the other eight and is
the wrong direction.  If the bug genuinely is in shared code, that is no
longer a debugging change — promote it to a unification change, fix it once
for everyone, and re-run `examples/baseline.sh` plus `tests/regress.sh` for
SPAsqr, SPAmix, SPACox, SPAGRM, SAGELD, WtCoxG, LEAF, SPAGxE and SPAGxEmix
before committing, with the per-column ledger the change requires.
