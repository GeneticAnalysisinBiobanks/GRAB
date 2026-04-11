# Genotype Factory

The `src/geno_factory/` directory contains all genotype I/O backends and
shared QC code.  Every statistical method in GRAB reads genotypes through
the abstract `GenoMeta` / `GenoCursor` interface defined in `geno_data.hpp`,
so adding or replacing a backend requires no changes to any analysis code.

## Architecture

```
┌─────────────────────────────────────────────────┐
│               marker engine                      │
│  (src/engine/marker.cpp)                         │
│  work-stealing thread pool, QC filters, output   │
└────────────────────┬────────────────────────────┘
                     │  GenoCursor::getGenotypes()
                     ▼
┌─────────────────────────────────────────────────┐
│             geno_data.hpp                        │
│  GenoMeta (metadata + cursor factory)            │
│  GenoCursor (per-thread genotype decoder)         │
│  GenoSpec / makeGenoData() / parseGenoIIDs()     │
└────────────────────┬────────────────────────────┘
                     │
       ┌─────────┬───┴────┬──────────┐
       ▼         ▼        ▼          ▼
   PlinkData  PgenData  VcfData  BgenData
   plink.hpp  pgen.hpp  vcf.hpp  bgen.hpp
   plink.cpp  pgen.cpp  vcf.cpp  bgen.cpp
       │         │        │          │
       └─────────┴───┬────┴──────────┘
                     │
                     ▼
              hwe.hpp / hwe.cpp
        HweExact()  statsFromCounts()
              GenoStats struct
```

### Admixture I/O

Ancestry-aware genotype reading (local-ancestry BED, MSP conversion)
is handled by `abed_io.hpp`, `abed_convert_txt.hpp`, and
`abed_convert_msp.hpp` in `src/localplus/`.

### Subject I/O

`src/subj_reader/` contains subject-level readers that are *not* genotype backends:

| File | Purpose |
|------|---------|
| `subject_data.hpp/.cpp` | Phenotype, covariate, and residual file I/O |
| `subject_filter.hpp/.cpp` | Subject inclusion/exclusion list parsing |
| `subject_set.hpp/.cpp` | Subject-set container (union/intersection masks) |
| `sparse_grm.hpp/.cpp` | Sparse GRM reader |

## File Inventory

| File | Lines | Role |
|------|------:|------|
| `geno_data.hpp` | ~110 | Abstract interface: `GenoMeta`, `GenoCursor`, `GenoSpec` |
| `geno_factory.cpp` | ~60 | `makeGenoData()` and `parseGenoIIDs()` dispatch |
| `hwe.hpp` | ~35 | `HweExact()`, `statsFromCounts()`, `GenoStats` declarations |
| `hwe.cpp` | ~130 | SNPHWE2 exact test implementation |
| `plink.hpp/.cpp` | ~620 | PLINK BED backend (plink2 SIMD decode + count) |
| `pgen.hpp/.cpp` | ~550 | PLINK2 PGEN backend (pgenlib wrapper) |
| `vcf.hpp/.cpp` | ~450 | VCF/BCF backend (htslib wrapper) |
| `bgen.hpp/.cpp` | ~400 | BGEN 1.2 backend (bgen reference library wrapper) |

## HWE Method

All backends use a single shared HWE exact test (`HweExact` in `hwe.cpp`),
implementing the SNPHWE2 algorithm:

> Wigginton JE, Cutler DJ, Abecasis GR (2005).
> "A Note on Exact Tests of Hardy-Weinberg Equilibrium."
> *Am J Hum Genet* 76:887–893.

This is the **plink2 `--hardy` default** method.  The chi-squared
approximation that was previously available has been removed — the exact
test is O(het_count) per marker with O(1) auxiliary memory and is
fast enough for all practical sample sizes.

### Why one method?

- The chi-squared approximation was inaccurate for rare variants
  (expected cell counts < 5).
- The previous code used a `bool exactHwe` switch that added complexity
  to every call path — the switch and all duplicate implementations
  have been removed.
- plink2 uses the exact test by default; there is no reason for GRAB to
  diverge.

## Third-Party Code Policy

The genotype factory wraps three external libraries:

| Library | Version | Source | Integration |
|---------|---------|--------|-------------|
| pgenlib | plink2 a.6.33 | `third_party/plink2-a.6.33/` | Compiled from source, headers in include path |
| htslib | 1.23.1 | `third_party/htslib-1.23.1/` | Compiled from source via its own Makefile |
| bgen | 1.2.0 | `third_party/bgen-1.2.0/` | Header-only reference implementation |

**Rules:**

1. Pin all third-party versions in `third_party/`.  Do not use system-installed
   libraries.
2. When upgrading, update the directory name to include the version number.
3. Do not modify third-party source files.  If a patch is needed, add a
   separate `.patch` file and document the reason.
4. Prefer wrapping over copying:  keep third-party APIs behind the
   `GenoCursor` interface so the rest of GRAB never calls them directly.

## Performance Design

### Bitmask-Based Subject Filtering

The `usedMask` (a `std::vector<uint64_t>` bitmask, one bit per sample in the
file) is passed to every backend at construction time.  Backends that support
it (PLINK, PGEN, BGEN) skip unused samples during decoding — the output
`Eigen::VectorXd` contains only used-subject genotypes in dense layout.

### Multithreading

The marker engine (`src/engine/marker.cpp`) uses a work-stealing thread pool:

- N worker threads pull chunks via `atomic<size_t>` counter
- 1 writer thread drains completed chunks in order
- Each worker owns a `GenoCursor` clone — no shared mutable state
- Per-thread buffers (genotype vector, missing indices, format buffer)
  are allocated once and reused

### PLINK BED Decoding (plink2 SIMD Primitives)

`plink.cpp` uses plink2's vectorized primitives from pgenlib for both
counting and decoding raw BED bytes:

- **All-used path:** `GenoarrCountFreqsUnsafe()` (carry-save-adder
  counting, word-at-a-time) + `GenoarrLookup16x8bx2()` (SIMD 16-entry
  lookup, 2 samples per 16-byte store) with a custom BED-code-to-double
  table (`kBedDoublePairs`).
- **Bitmask path:** `GenoarrCountSubsetFreqs2()` (popcount-based subset
  counting) + scalar bit-scan scatter for decode.
- pgenlib handles SIMD tier selection automatically: SSE2 (x86-64
  baseline), AVX2 (when `-mavx2` is set), or ARM via SIMDe.
- Raw BED bytes are copied to a word-aligned scratch buffer
  (`m_alignedBed`) before calling plink2 functions.
