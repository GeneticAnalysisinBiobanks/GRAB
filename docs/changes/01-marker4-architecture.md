# Marker4 Architecture Changes

## Scope
This file summarizes the new `GRAB.Marker4` + `Main4.cpp` architecture and API behavior.

## R-side entry point
- Added `R/GRAB_Marker4.R`.
- New user-facing function: `GRAB.Marker4(...)`.
- Supports method classes:
  - `POLMM_NULL_Model`
  - `SPACox_NULL_Model`
  - `SPAmix_NULL_Model`
  - `SPAmixPlus_NULL_Model`
  - `SPAGRM_NULL_Model`
  - `SAGELD_NULL_Model`
  - `WtCoxG_NULL_Model`
  - `SPAsqr_NULL_Model`
  - `LEAF_NULL_Model`

## New controls and behavior
- `nthreads` argument controls C++ worker count.
- `nthreads = NULL` uses `data.table::getDTthreads()`.
- `overwrite` argument controls whether existing output is removed.
- `control$omp_num_threads` is rejected for `GRAB.Marker4`; thread control is centralized in `nthreads`.
- Added runtime logs including genotype type and active thread count.

## Chunk execution model
- R constructs `genoIndexList` by chromosome and `nMarkersEachChunk`.
- R passes chunks in a single call to C++ (`mainMarkerChunksInCPP4`) instead of looping in R.
- C++ performs chunk scheduling, marker testing, and output writing.

## C++ implementation (`src/Main4.cpp`)
- Added `mainMarkerChunksInCPP4(...)` with merged marker/QC params and optional `extraParams`.
- Uses dynamic chunk scheduling with `std::atomic<size_t> nextChunk`.
- Uses one writer thread and multiple worker threads.
- Preserves deterministic output order by chunk index.
- Uses per-thread genotype readers (`PLINK::PlinkClass` / `BGEN::BgenClass`) to avoid reader contention.

## Method-specific output compatibility
`Main4.cpp` keeps existing output schema semantics by method:
- POLMM
- SPACox
- SPAmix / SPAmixPlus
- SPAGRM
- SAGELD/GALLOP
- WtCoxG
- SPAsqr
- LEAF

## Export wiring
- Added C++ export/wrapper for `mainMarkerChunksInCPP4`.
- Updated registration and R wrapper signature in generated `RcppExports` files.

## Shared object access across translation units
- Exposed globals from `Main.cpp` through `Main.h` (`extern`) so `Main4.cpp` can reuse initialized objects.
- This was required because file-local `static` globals in `Main.cpp` are not visible from `Main4.cpp`.
