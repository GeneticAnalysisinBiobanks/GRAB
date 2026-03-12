# Validation Summary

## Test script
- Executed: `Rscript z.test.R`
- Script flow:
  1. Fit `SPAsqr` null model.
  2. Run baseline `GRAB.Marker`.
  3. Run `GRAB.Marker4` with `nthreads = 1`.
  4. Run `GRAB.Marker4` with `nthreads = 10`.
  5. Compare result tables using `all.equal`.

## Outcomes
- Build and package load succeeded.
- `GRAB.Marker4` completed for both 1-thread and 10-thread runs.
- Equality checks in script returned `TRUE` for comparisons to baseline in current test.
- No stack imbalance warnings in the latest post-fix run.
- No segmentation fault in latest post-fix run.

## Notes
- Compiler warnings from Boost templates were present during build, but they were not runtime-fatal and did not block analysis completion.
- Current behavior indicates that the SPAsqr multithread path is functionally stable for this regression case.

## Suggested follow-up validation
- Add a larger PLINK/BGEN benchmark to measure speedup scaling (`nthreads = 1, 2, 4, 8, 16`).
- Add CI-style regression for `SPAsqr` threaded output equivalence.
- Add stress test with many chunks and mixed missingness patterns.
