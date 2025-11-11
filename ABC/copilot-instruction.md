# GitHub Copilot Instructions for GRAB

## Overview

GRAB is a high-performance R package for genome-wide association studies (GWAS) in biobank-scale data. It implements multiple statistical methods as a hybrid R/C++ architecture with method-specific optimizations for different phenotypes.

## Architecture & Core Concepts

### Two-Step GWAS Workflow
All methods follow this pattern:
1. **Step 1**: `GRAB.NullModel()` - Fit null model, estimate parameters, handle relatedness
2. **Step 2**: `GRAB.Marker()` / `GRAB.Region()` - Conduct association testing using null model object

### Statistical Methods (src/[METHOD].cpp + R/[METHOD].R)
- **POLMM**: Ordinal categorical traits with relatedness (`"ordinal"`)
- **SPACox**: Time-to-event for unrelated samples (`"time-to-event"`)  
- **SPAmix**: Time-to-event for admixed populations (`"time-to-event"` + PCs)
- **SPAGRM**: Complex traits with relatedness (`"Residual"`)
- **WtCoxG**: Case-ascertained survival analysis (`"time-to-event"`)
- **SAGELD**: Longitudinal/complex traits (`"Residual"`)

### File Structure Patterns
```
R/
├── GRAB_Null_Model.R      # Main null model fitting entry point
├── GRAB_Marker.R          # Single-variant testing entry point  
├── GRAB_Region.R          # Gene-based testing entry point
├── [METHOD].R             # Method-specific R implementations
├── control.R              # Parameter validation & defaults
├── Geno.R                 # Genotype file handling (PLINK/BGEN)
└── Util.R                 # Shared utilities & QC functions

src/
├── Main.cpp               # Primary R-C++ interface & global objects
├── [METHOD].cpp/.h        # Method-specific C++ implementations  
├── PLINK.cpp/.h           # PLINK format handlers
├── BGEN.cpp/.h            # BGEN format handlers
└── RcppExports.cpp        # Auto-generated R-C++ bindings
```

## Critical Development Patterns

### Method Organization
Each method implements this interface:
```r
# R side (in R/[METHOD].R)
checkControl.NullModel.[METHOD](control, traitType)
checkControl.Marker.[METHOD](control) 
checkControl.Region.[METHOD](control)
fitNullModel.[METHOD](response, designMat, subjData, control, ...)
setMarker.[METHOD](objNull, control)
mainMarker.[METHOD](genoType, genoIndex, outputColumns, objNull)
```

```cpp
// C++ side (in src/[METHOD].cpp) 
class [METHOD]Class {
  void getMarkerPval(arma::vec& GVec, ...);  // Core association test
  arma::vec getpvalVec();                    // Return p-values
};
```

### Global State Management (src/Main.cpp)
- **Global objects**: `ptr_g[METHOD]obj` store method-specific state
- **Setup functions**: `set[METHOD]objInCPP()` initialize before analysis
- **Config globals**: `g_impute_method`, `g_missingRate_cutoff`, etc.

### Control Parameter Philosophy
Replace `eval(parse())` patterns with explicit if-else chains:
```r
# Replace this:
textToParse <- paste0("func.", method, "(args)")
eval(parse(text = textToParse))

# With this:
if (NullModelClass == "POLMM_NULL_Model") {
  result <- func.POLMM(args)
} else if (NullModelClass == "SPACox_NULL_Model") {
  result <- func.SPACox(args)
} # ... etc
```

### Genotype File Handling
- **Unified interface**: `Unified_getOneMarker()` abstracts PLINK/BGEN
- **Memory management**: Process in chunks via `nMarkersEachChunk`
- **Subject filtering**: `subjData` parameter controls sample inclusion
- **Missing data**: Three strategies: `"mean"`, `"minor"`, `"drop"`

## Build & Development Workflow

### Package Build
```bash
# Standard R package commands
R CMD build .
R CMD check --as-cran GRAB_*.tar.gz
R CMD INSTALL GRAB_*.tar.gz
```

### C++ Compilation Dependencies
- **Makevars**: Links against LAPACK, BLAS, zlib with Armadillo 64-bit indexing
- **External libs**: RcppArmadillo, BH (Boost), data.table
- **Header structure**: Each method has `.h` declaring class interface

### Example Testing Pattern
All methods include roxygen examples following this pattern:
```r
#' @examples
#' # Step 1: Fit null model
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' obj.METHOD <- GRAB.NullModel(
#'   formula, data = PhenoData, subjData = IID,
#'   method = "METHOD", traitType = "...", ...
#' )
#'
#' # Step 2: Association testing
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "Results.txt")
#' GRAB.Marker(obj.METHOD, GenoFile = GenoFile, OutputFile = OutputFile)
```

## Integration Points & Dependencies

### R Dependencies
- **Core analysis**: data.table, survival, SKAT 
- **C++ integration**: Rcpp, RcppArmadillo
- **Optional**: Matrix, CompQuadForm for specific methods

### External File Formats
- **PLINK**: `.bed/.bim/.fam` binary format (most common)
- **BGEN**: `.bgen` + `.bgi` index for UK Biobank data
- **Sparse GRM**: Custom tab-delimited format (`ID1 ID2 Value`)

### Output File Management  
- **Chunked processing**: Large analyses split across multiple chunks
- **Resume capability**: `.index` files track progress for restarts
- **Multiple outputs**: Main results + `.markerInfo` + `.otherMarkerInfo`

## Performance Considerations

### Memory & Threading
- **OpenMP**: `omp_num_threads` controls C++ parallelization
- **Memory chunks**: `memoryChunk` limits RAM usage during genotype reading
- **QC thresholds**: `missing_cutoff`, `min_maf_marker`, `min_mac_marker` filter variants

### Method-Specific Optimizations
- **POLMM**: Dense/Sparse GRM options, LOCO analysis for relatedness
- **SPAmix**: Requires PC columns for population structure
- **WtCoxG**: Two-step process with batch effect testing
- **Region analysis**: SKAT-O, burden tests, ultra-rare variant collapsing

Focus on these patterns when working with method implementations, parameter validation, or extending functionality. The codebase prioritizes type safety, performance, and robust statistical computation over generic programming patterns.
