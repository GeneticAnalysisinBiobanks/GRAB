# Release Notes

## v0.4.0 (2026-03-25)

  - Null object data is transferred from R to C++ in a single handoff.
  - All methods follow consistent low-level optimizations for speed and memory.
  - Computation is parallelized by chunk for maximum scalability on multicore systems.

## v0.3.0 (2026-02-06)

- Added new methods:
  - SAGELD, developed by He Xu
  - LEAF, developed by Ying Li
  - Added SPAsqr, developed by Dr. Heng
  - SPAmixPlus and SPAmixLocalPlus, developed by Yuzhuo Ma

- Changes to SPAsqr:
  - Added quantile regression
  - Changed outlier cutoff from 0.7 to 0.8
  - Added SPAsqr.Step1a, which builds GRM cache for all families
  - Added SPAsqr.Step1b, which fits null model using pre-computed GRM cache

- Changes to existing methods:
  - Added verbose options (default FALSE) to functions in Geno.R
  - WtCoxG now returns Z-scores and score statistics in addition to p-values

## v0.2.4 (2025-12-01)

- Reorganized R code; helper functions executed only once are now inlined.
- Revised roxygen2 documentation to comprehensively list all available parameters.
- Refactored `GRAB.NullModel()`, `GRAB.Marker()`, and `Main.cpp`.
- Reintegrated WtCoxG with improved structure.

## v0.2.3 (2025-08-25)

- Reformatted and significantly improved logging output.
- Refined function help documentation and examples.
- Move the `TestforBatchEffect()` call from `GRAB.NullModel()` to `fitNullModel.WtCoxG()`.
- Update all examples so that they never write to the installation path.

## v0.2.2 (2025-08-01)

The following notes summarize the main changes from v0.2.0 to v0.2.2.

Resolved all issues raised by the CRAN manual review to meet publication requirements, including:

- Updated the `Description` field in `DESCRIPTION` to a single-paragraph summary with relevant references.
- Added a `@return` section to the documentation of all exported functions where absent.
- Replaced all `print()` and `cat()` statements with `message()` to allow users to suppress output if desired.
- Removed code that set specific random seeds inside functions.
- Redirect example output files to the R session's `tempdir()` to avoid side effects on the user's file system.

Increased the number of markers (which was decreased from 10k to 1k to minimize package size and running time) in the example dataset from 1k to 1.1k to ensure the number of markers exceeds the number of subjects.

Moved inst/docker back to the project root and excluded it from the package build.

## v0.2.0 (2025-07-17)

Add a new method, WtCoxG, and bump version to 0.2.0.

We reorganized the `Depends`, `Imports`, `Suggests` and `LinkingTo` fields in the `DESCRIPTION` file for clarity and correctness. Additionally, we consolidated `Makevars` and `Makevars.win` into a single `Makevars` file for simplified build configuration.

In previous versions, GRAB automatically imported all functions and objects from its dependent packages and exported all package objects to the user's workspace. We have updated the NAMESPACE file so that GRAB now only imports specific functions as needed from dependencies and only exports user-facing functions. This results in a cleaner namespace and avoids cluttering the user's environment with unnecessary objects.

We no longer include a copy of GCTA program within the package. Instead, users have to specify the path of the GCTA executable on their system using `gcta64File` parameter when calling `getTempFilesFullGRM()`. We have successfully tested `getTempFilesFullGRM()` with GCTA v1.94.4 on Windows. Despite previous documentation indicating Linux-only support, this function works on Windows systems as well.

We have tested and debugged all examples included in the function help files, and verified that they execute correctly. We have minimized the package size and example running time. Additionally, we have made all necessary changes to ensure the package complies with CRAN policies and requirements, passing `R CMD check --as-cran` with no `ERROR` or `WARNING`, and as few `NOTE` as possible.

The old version 0.1.2 can be found in branch [release/v0.1.2](https://github.com/GeneticAnalysisinBiobanks/GRAB/tree/release/v0.1.2). This branch will serve as an archive and will no longer be actively maintained.
