# Thread-Safety and SPAsqr Refactor

## Problem observed
- `GRAB.Marker4` with `SPAsqr` and high thread counts produced runtime instability:
  - stack/SEXP-related errors (`SET_STRING_ELT` / stack imbalance)
  - segmentation fault in earlier runs

## Root causes addressed
- R API usage from worker-thread setup paths.
- Risky ownership/copy patterns in SPAsqr-related objects under multithread fan-out.
- Remaining Rcpp-heavy structures in SPAGRM hot paths used by SPAsqr marker testing.

## Fixes applied

### 1) SPAGRM hot-path storage conversion
Files:
- `src/SPAGRM.h`
- `src/SPAGRM.cpp`

Changes:
- Replaced Rcpp list members with native C++ containers:
  - `m_TwoSubj_resid_list`
  - `m_TwoSubj_rho_list`
  - `m_ThreeSubj_standS_list`
  - `m_ThreeSubj_CLT_list`
- Constructor now converts incoming R lists once.
- Marker-time calculations (`MGF_cpp`, root finding, SPA probability path) use native containers.

### 2) SPAsqr object ownership update
File:
- `src/SPAsqr.h`

Changes:
- Replaced `std::vector<SPAGRM::SPAGRMClass*>` with `std::vector<SPAGRM::SPAGRMClass>`.
- Removed manual destructor and pointer deletes.
- Uses value semantics (`emplace_back`) and direct access (`.`) to avoid pointer lifetime/copy hazards.

### 3) Reader setup made worker-safe
Files:
- `src/PLINK.cpp`
- `src/BGEN.cpp`

Changes:
- Replaced `Rcpp::match` sample index mapping with pure C++ hash-map matching.
- Removed worker-thread dependence on R API during reader construction.

### 4) SPAsqr multithreading re-enabled
File:
- `R/GRAB_Marker4.R`

Changes:
- Removed `SPAsqr_NULL_Model` from temporary forced-single-thread fallback.
- Kept fallback only for methods still marked unsafe (`WtCoxG`, `LEAF`).

## Net effect
- SPAsqr in `GRAB.Marker4` now runs stably at `nthreads > 1` in current regression checks.
- Worker initialization and SPAsqr marker paths are now aligned with thread-safe usage expectations.
