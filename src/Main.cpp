/**
 * @file Main.cpp
 * @brief Main interface between C++ and R for the GRAB package
 *
 * This file contains the core C++ functions that are exposed to R via Rcpp.
 * It serves as the main entry point for genetic association analysis methods
 * including POLMM, POLMM-GENE, SPACox, SPAmix, SPAGRM, SAGELD and WtCoxG.
 *
 * Functions (one‑line summaries):
 *   setSparseGRMInCPP            : Load sparse GRM from R list (triplet form) into Armadillo.
 *   setDenseGRMInCPP             : Construct / reset dense GRM object with QC thresholds.
 *   getDenseGRMInCPP             : Multiply dense GRM by a vector.
 *   setMarker_GlobalVarsInCPP    : Set global QC + threading + grouping parameters for marker tests.
 *   setRegion_GlobalVarsInCPP    : Set global parameters for region / rare variant analyses.
 *   updateGroupInfo              : Compute per‑group sample counts, alt counts and frequencies.
 *   getLabelInfo                 : Compute per‑label MAC and MAF matrices (case/control or multi‑group).
 *   Unified_getOneMarker         : Fetch one marker's genotypes + metadata from PLINK/BGEN.
 *   Unified_getMarkerPval        : Compute single‑marker p‑value plus HWE p‑value.
 *   Unified_getRegionPVec        : Region / aggregate test statistics (normal + SPA forms).
 *   mainMarkerInCPP              : High‑level GWAS loop over selected markers (multi‑phenotype aware).
 *   mainRegionURVInCPP           : Collapse ultra‑rare variants and return burden statistics.
 *   mainRegionInCPP              : Comprehensive gene/region analysis (single, burden, URV, VC tests).
 *   getGenoInfoInCPP             : Quick genotype summary (freq/missing) without full matrix load.
 *   getGenoInCPP                 : Load genotype matrix (adaptive marker count) with imputation.
 *   getGenoInCPP_fixedNumber     : Load fixed number of markers with missing/MAF filtering.
 *   getSpGenoInCPP               : Load sparse genotype matrix representation.
 *   setPLINKobjInCPP             : Initialize PLINK reader object for downstream access.
 *   setBGENobjInCPP              : Initialize BGEN reader object with sample filtering.
 *   closeGenoInputInCPP          : Release genotype reader resources (PLINK/BGEN).
 *   setPOLMMobjInCPP             : Initialize POLMM object with pre‑computed matrices.
 *   setPOLMMobjInCPP_NULL        : Fit POLMM null model and prepare C++ state.
 *   setSPAGRMobjInCPP            : Initialize SPAGRM object with residual / GRM decomposition.
 *   setSAGELDobjInCPP            : Initialize SAGELD/GALLOP object with residual partitions.
 *   setSPAmixobjInCPP            : Initialize SPAmix object (multi‑phenotype residual + PCs).
 *   setSPACoxobjInCPP            : Initialize SPACox object (CGF grid + residual info).
 *   setWtCoxGobjInCPP            : Initialize weighted Cox (WtCoxG) object with objNull$mergeGenoInfo.
 */

// Rcpp dependencies for statistical computing and linear algebra
// [[Rcpp::depends(BH)]]          // Boost headers for advanced algorithms
// [[Rcpp::depends(RcppArmadillo)]]  // Armadillo C++ linear algebra library
#include <RcppArmadillo.h>
#include <boost/math/distributions/beta.hpp>  // Beta distribution functions

// Standard C++ libraries
#include <thread>   // std::this_thread::sleep_for - for thread management
#include <chrono>   // std::chrono::seconds - for time-based operations
#include <cstdio>   // std::remove - for file operations

// Include headers for various components of the GRAB package
#include "PLINK.h"     // PLINK format genotype file handler
#include "BGEN.h"      // BGEN format genotype file handler
#include "POLMM.h"     // Polytomous Logistic Mixed Model analysis
#include "UTIL.h"      // Utility functions for statistical computations
#include "SPACox.h"    // Score test with Pooled Approximate likelihood (Cox regression)
#include "DenseGRM.h"  // Dense Genetic Relationship Matrix operations
#include "SPAmix.h"    // Mixed effect models with SPA correction
#include "SPAGRM.h"    // SPA test with GRM (Genetic Relationship Matrix)
#include "SAGELD.h"    // Scalable and Accurate Genomic analysis for Extreme Large Data
#include "WtCoxG.h"    // Weighted Cox regression for genetic data


//==============================================================================
// SECTION 1: OBJECT DECLARATIONS
//==============================================================================

// Global object pointers for different genotype file formats
static PLINK::PlinkClass* ptr_gPLINKobj = nullptr;   // PLINK format
static BGEN::BgenClass* ptr_gBGENobj = nullptr;      // BGEN format
// static VCF::VcfClass* ptr_gVCFobj = nullptr;       // VCF format (disabled for now)

// Global object for dense genetic relationship matrix operations
static DenseGRM::DenseGRMClass* ptr_gDenseGRMobj = nullptr;

// Global object pointers for different statistical analysis methods
static POLMM::POLMMClass* ptr_gPOLMMobj = nullptr;
static SPACox::SPACoxClass* ptr_gSPACoxobj = nullptr;
static SPAmix::SPAmixClass* ptr_gSPAmixobj = nullptr;
static SPAGRM::SPAGRMClass* ptr_gSPAGRMobj = nullptr;
static SAGELD::SAGELDClass* ptr_gSAGELDobj = nullptr;
static WtCoxG::WtCoxGClass* ptr_gWtCoxGobj = nullptr;

// Global configuration variables for genetic analysis
static std::string g_impute_method;          // Imputation method: "mean", "minor", or "drop"
static double g_missingRate_cutoff;          // Maximum allowed missing rate for markers
static unsigned int g_omp_num_threads;       // Number of OpenMP threads for parallel processing

// Quality control thresholds for marker-level analysis
static double g_marker_minMAF_cutoff;         // Minimum Minor Allele Frequency for single markers
static double g_marker_minMAC_cutoff;         // Minimum Minor Allele Count for single markers

// Quality control thresholds for region-based analysis
static double g_region_minMAC_cutoff;  // Min MAC for rare variants aggregation (like SAIGE-GENE+)
static double g_region_maxMAF_cutoff;         // Maximum MAF to consider variants as "rare"
static unsigned int g_region_maxMarkers_cutoff;  // Max markers per chunk (memory management)
static arma::vec g_region_weight_beta;        // Beta parameters for variant weighting
static arma::vec g_region_max_maf_vec;        // MAF thresholds for different variant categories

// Group-specific variables (used only for POLMM analysis)
static arma::uvec g_group;                    // Group assignment for each individual
static bool g_ifOutGroup;                     // Whether to output group-specific statistics
static unsigned int g_nGroup;                 // Total number of groups

// Sparse Genetic Relationship Matrix (GRM) for kinship correction
static arma::sp_mat g_SparseGRM;

// Performance monitoring variables
static arma::vec g_compTime1(2, arma::fill::zeros);  // Timing for Unified_getOneMarker function
static arma::vec g_compTime2(2, arma::fill::zeros);  // Timing for Unified_getRegionPVec function
static arma::vec g_compTime3(2, arma::fill::zeros);  // Additional timing measurements


//==============================================================================
// SECTION 2: GLOBAL CONFIGURATION FUNCTIONS
// Functions to set up analysis parameters and data structures
//==============================================================================

// Converts an R sparse matrix representation into a C++ Armadillo sparse matrix
// for efficient kinship-based corrections in genetic association analysis.
// [[Rcpp::export]]
void setSparseGRMInCPP(
  Rcpp::List t_KinMatListR  // R list containing sparse kinship matrix data
) {
  // Extract sparse matrix components from R list
  arma::umat locations = t_KinMatListR["locations"];  // Row/column indices of non-zero elements
  arma::vec values = t_KinMatListR["values"];         // Values at those locations
  int n = t_KinMatListR["nSubj"];                     // Number of subjects (matrix dimension)

  // Construct sparse matrix from triplet format (locations, values, dimensions)
  arma::sp_mat KinMat(
    locations,  // Row/column indices of non-zero elements
    values,     // Values at those locations
    n,          // Number of rows
    n           // Number of columns
  );
  g_SparseGRM = KinMat;  // Store in global variable
}

// Creates a dense genetic relationship matrix object that can efficiently
// compute kinship corrections for large-scale genetic association analysis.
// [[Rcpp::export]]
void setDenseGRMInCPP(
  double t_memoryChunk,     // Memory allocation size for dense GRM computation (in GB)
  double t_minMafGRM,       // Minimum MAF threshold for variants used in GRM construction
  double t_maxMissingGRM    // Maximum missing rate allowed for variants in GRM
) {
  // Clean up existing dense GRM object if it exists
  if (ptr_gDenseGRMobj)
    delete ptr_gDenseGRMobj;

  // Create new dense GRM object with specified parameters
  ptr_gDenseGRMobj = new DenseGRM::DenseGRMClass(
    ptr_gPLINKobj,     // PLINK genotype file reader object
    t_memoryChunk,     // Memory allocation size (in GB)
    t_minMafGRM,       // Minimum MAF threshold for GRM construction
    t_maxMissingGRM    // Maximum missing rate for GRM variants
  );
}

// Performs efficient matrix-vector multiplication K*b where K is the dense GRM
// and b is the input vector. Used for Leave-One-Chromosome-Out (LOCO) analysis.
// [[Rcpp::export]]
arma::vec getDenseGRMInCPP(
  arma::vec t_bVec,         // Input vector to be corrected (typically residuals or phenotypes)
  std::string t_excludeChr, // Chromosome to exclude from kinship computation (for LOCO analysis)
  int t_grainSize           // Parallel processing granularity parameter
) {
  arma::vec yVec = DenseGRM::getKinbVec(
    t_bVec,           // Input vector to be corrected
    ptr_gDenseGRMobj, // Dense GRM object pointer
    t_excludeChr,     // Chromosome to exclude for LOCO analysis
    t_grainSize       // Parallel processing granularity
  );
  return yVec;
}

// Sets up quality control parameters and computational settings for single-marker
// genome-wide association studies.
// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(
  std::string t_impute_method,    // Method for handling missing genotypes ("mean", "minor", "drop")
  double t_missing_cutoff,        // Maximum allowed missing rate for markers
  double t_min_maf_marker,        // Minimum Minor Allele Frequency threshold
  double t_min_mac_marker,        // Minimum Minor Allele Count threshold
  unsigned int t_omp_num_threads  // Number of OpenMP threads for parallel computation
) {
  // Set global configuration variables
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_omp_num_threads = t_omp_num_threads;
}

// Configures parameters for gene-based rare variant association analysis,
// following approaches similar to SAIGE-GENE+ for burden and variance-component tests.
// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(
  std::string t_impute_method,      // Method for handling missing genotypes
  double t_missing_cutoff,          // Maximum allowed missing rate for variants
  double t_max_maf_region,          // Maximum MAF to consider variants as "rare"
  double t_min_mac_region,          // Minimum MAC threshold for rare variant aggregation
  unsigned int t_max_markers_region, // Maximum number of markers per analysis chunk
  unsigned int t_omp_num_threads,   // Number of OpenMP threads for parallel computation
  arma::vec t_region_weight_beta,   // Beta parameters for variant weighting schemes
  arma::vec t_region_max_maf_vec    // MAF thresholds for different variant categories
) {
  // Set basic quality control parameters
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_omp_num_threads = t_omp_num_threads;

  // Set region-specific analysis parameters
  g_region_minMAC_cutoff = t_min_mac_region;
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_region_maxMarkers_cutoff = t_max_markers_region;

  // Set variant weighting parameters
  g_region_weight_beta = t_region_weight_beta;
  g_region_max_maf_vec = t_region_max_maf_vec;
}

//==============================================================================
// SECTION 3: CORE UTILITY FUNCTIONS
// Helper functions used by main analysis routines
//==============================================================================

// Computes group-stratified allele frequency and count statistics, accounting
// for missing data patterns. Used in population stratification analysis.
void updateGroupInfo(
  arma::vec t_GVec,                        // Genotype vector for current marker
  std::vector<uint32_t> t_indexForMissing, // Indices of individuals with missing genotypes
  arma::vec &nSamplesInGroupVec,           // Output: number of non-missing samples per group
  arma::vec &AltCountsInGroupVec,          // Output: alternate allele counts per group
  arma::vec &AltFreqInGroupVec             // Output: alternate allele frequencies per group
) {
  unsigned int n1 = t_GVec.size();
  nSamplesInGroupVec.zeros();  // Initialize group sample counts
  AltCountsInGroupVec.zeros();  // Initialize group allele counts

  // Handle edge case: if no missing data, add sentinel value
  if (t_indexForMissing.size() == 0)
    t_indexForMissing.push_back(n1);

  // Iterate through all individuals, skipping those with missing genotypes
  unsigned int i1 = 0;  // Index into missing data vector
  for (unsigned int i = 0; i < n1; i++) {
    // Check if current individual has missing genotype
    if (i == t_indexForMissing.at(i1)) {
      // Move to next missing index if available
      if (i1 < t_indexForMissing.size() - 1)
        i1++;
    } else {
      // Individual has non-missing genotype - update group statistics
      unsigned int grp = g_group.at(i);  // Get group assignment for this individual

      nSamplesInGroupVec.at(grp) += 1;              // Increment sample count for this group
      AltCountsInGroupVec.at(grp) += t_GVec.at(i);  // Add genotype to allele count
    }
  }

  // Calculate allele frequencies from counts (divide by 2 for diploid genotypes)
  AltFreqInGroupVec = AltCountsInGroupVec / nSamplesInGroupVec / 2;
}

// Computes label-stratified allele statistics for case-control or
// multi-group analysis designs.
void getLabelInfo(
  arma::vec t_GVec,                    // Genotype vector after imputation
  std::vector<unsigned int> t_labelVec, // Label assignment for each individual
  unsigned int t_rowIndex,             // Row index in output matrices
  arma::mat &t_MACLabelMat,            // Output matrix for Minor Allele Counts by label
  arma::mat &t_MAFLabelMat             // Output matrix for Minor Allele Frequencies by label
) {
  unsigned int n = t_labelVec.size();
  unsigned int nLabel = t_MACLabelMat.n_cols;

  arma::vec MACLabelVec(nLabel, arma::fill::zeros);
  arma::vec counttLabelVec(nLabel, arma::fill::zeros);

  // Count samples and alleles per label
  for (unsigned int i = 0; i < n; i++) {
    counttLabelVec.at(t_labelVec.at(i) - 1) += 1;            // Sample count per label
    MACLabelVec.at(t_labelVec.at(i) - 1) += t_GVec.at(i);    // Allele count per label
  }

  // Calculate allele frequencies (divide by 2 for diploid)
  arma::vec MAFLabelVec = MACLabelVec / (counttLabelVec * 2);

  // Store results in output matrices
  t_MACLabelMat.row(t_rowIndex) = MACLabelVec.t();
  t_MAFLabelMat.row(t_rowIndex) = MAFLabelVec.t();
}

//==============================================================================
// SECTION 4: CORE GENOTYPE AND STATISTICAL ANALYSIS FUNCTIONS
// Main computational functions for genetic association analysis
//==============================================================================

// Unified interface for extracting genotype data from different file formats.
// Handles missing data tracking and provides comprehensive marker information.
arma::vec Unified_getOneMarker(
  std::string t_genoType,                    // Genotype file format ("PLINK", "BGEN")
  uint64_t t_gIndex,                         // Marker index in genotype file
  std::string &t_ref,                        // Reference allele (output)
  std::string &t_alt,                        // Alternate allele (output)
  std::string &t_marker,                     // Marker ID (output)
  uint32_t &t_pd,                            // Physical position (output)
  std::string &t_chr,                        // Chromosome (output)
  double &t_altFreq,                         // Alternate allele frequency (output)
  double &t_altCounts,                       // Alternate allele count (output)
  double &t_missingRate,                     // Missing genotype rate (output)
  double &t_imputeInfo,                      // Imputation quality score (output)
  bool t_isOutputIndexForMissing,            // Whether to output missing data indices
  std::vector<uint32_t> &t_indexForMissing,  // Indices of missing genotypes (output)
  bool t_isOnlyOutputNonZero,                // Whether to output only non-zero genotypes
  std::vector<uint32_t> &t_indexForNonZero   // Indices of non-zero genotypes (output)
) {
  arma::vec GVec;

  if (t_genoType == "PLINK") {
    GVec = ptr_gPLINKobj->getOneMarker(
      t_gIndex,                        // Marker index in genotype file
      t_ref,                           // Reference allele (output)
      t_alt,                           // Alternate allele (output)
      t_marker,                        // Marker ID (output)
      t_pd,                            // Physical position (output)
      t_chr,                           // Chromosome (output)
      t_altFreq,                       // Alternate allele frequency (output)
      t_altCounts,                     // Alternate allele count (output)
      t_missingRate,                   // Missing genotype rate (output)
      t_imputeInfo,                    // Imputation quality score (output)
      t_isOutputIndexForMissing,       // Whether to output missing data indices
      t_indexForMissing,               // Indices of missing genotypes (output)
      t_isOnlyOutputNonZero,           // Whether to output only non-zero genotypes
      t_indexForNonZero,               // Indices of non-zero genotypes (output)
      true                             // Is dosage read flag
    );
  }

  if (t_genoType == "BGEN") {
    bool t_isBoolRead = false;
    GVec = ptr_gBGENobj->getOneMarker(
      t_gIndex,                        // Marker index in genotype file
      t_ref,                           // Reference allele (output)
      t_alt,                           // Alternate allele (output)
      t_marker,                        // Marker ID (output)
      t_pd,                            // Physical position (output)
      t_chr,                           // Chromosome (output)
      t_altFreq,                       // Alternate allele frequency (output)
      t_altCounts,                     // Alternate allele count (output)
      t_missingRate,                   // Missing genotype rate (output)
      t_imputeInfo,                    // Imputation quality score (output)
      t_isOutputIndexForMissing,       // Whether to output missing data indices
      t_indexForMissing,               // Indices of missing genotypes (output)
      t_isOnlyOutputNonZero,           // Whether to output only non-zero genotypes
      t_indexForNonZero,               // Indices of non-zero genotypes (output)
      t_isBoolRead                     // Boolean read flag
    );
  }

  return GVec;
}

// Performs region-based association testing for rare variant analysis.
// Computes both normal approximation and SPA-corrected p-values.
void Unified_getRegionPVec(
  std::string t_method,    // Statistical method
  arma::vec t_GVec,        // Aggregated genotype vector for region
  double &t_Stat,          // Score statistic (output)
  double &t_Beta,          // Effect size estimate (output)
  double &t_seBeta,        // Standard error of effect size (output)
  double &t_pval0,         // Normal approximation p-value (output)
  double &t_pval1,         // SPA-corrected p-value (output)
  arma::vec &t_P1Vec,      // First component vector for variance calculation (output)
  arma::vec &t_P2Vec       // Second component vector for variance calculation (output)
) {
  if (t_method == "POLMM") {
    ptr_gPOLMMobj->getRegionPVec(
      t_GVec,   // Aggregated genotype vector for region
      t_Stat,   // Score statistic (output)
      t_Beta,   // Effect size estimate (output)
      t_seBeta, // Standard error of effect size (output)
      t_pval0,  // Normal approximation p-value (output)
      t_pval1,  // SPA-corrected p-value (output)
      t_P1Vec,  // First component vector for variance calculation (output)
      t_P2Vec   // Second component vector for variance calculation (output)
    );
  } else {
    Rcpp::stop("Internal errer in Unified_getRegionPVec: Method " + t_method);
  }
}

//==============================================================================
// SECTION 5: MAIN ANALYSIS FUNCTIONS
// High-level entry points for marker and region-based association analysis
//==============================================================================

// [[Rcpp::export]]
Rcpp::List mainMarkerInCPP(
  std::string t_method,                 // Statistical analysis method
  std::string t_genoType,               // Genotype file format
  std::vector<uint64_t> t_genoIndex     // Marker indices to analyze
) {
  int q = t_genoIndex.size(); // Number of markers to analyze

  // Initialize output vectors for marker information and statistics
  std::vector<std::string> markerVec(q);    // Marker IDs (e.g., rs12345)
  std::vector<std::string> infoVec(q);      // Marker info (CHR:POS:REF:ALT format)
  std::vector<double> altFreqVec(q);        // Alternate allele frequencies
  std::vector<double> altCountsVec(q);      // Alternate allele counts
  std::vector<double> missingRateVec(q);    // Missing genotype rates
  std::vector<double> hwepvalVec(q, arma::datum::nan);  // Hardy-Weinberg p-values

  // Determine number of phenotypes based on analysis method
  int Npheno = 1;  // Default: single phenotype
  if (t_method == "SPAmix")
    Npheno = ptr_gSPAmixobj->getNpheno();  // Mixed effects: multiple phenotypes
  if (t_method == "SAGELD")
    Npheno = 2;  // SAGELD: typically binary trait analysis
  if (t_method == "WtCoxG")
    Npheno = 2;  // WtCoxG: returns two p-values (with and without external reference)

  // Initialize test result vectors (sized for multiple phenotypes)
  std::vector<double> pvalVec(q * Npheno, arma::datum::nan);     // P-values
  std::vector<double> zScoreVec(q * Npheno, arma::datum::nan);   // Z-scores
  std::vector<double> BetaVec(q * Npheno, arma::datum::nan);     // Effect sizes (beta)
  std::vector<double> seBetaVec(q * Npheno, arma::datum::nan);   // Standard errors

  // Initialize group-stratified analysis matrices (if requested)
  arma::mat nSamplesInGroup;   // Sample counts per group
  arma::mat AltCountsInGroup;  // Allele counts per group
  arma::mat AltFreqInGroup;    // Allele frequencies per group

  if (g_ifOutGroup) {
    nSamplesInGroup.resize(
      q,        // Number of markers
      g_nGroup  // Number of groups
    );
    AltCountsInGroup.resize(
      q,        // Number of markers
      g_nGroup  // Number of groups
    );
    AltFreqInGroup.resize(
      q,        // Number of markers
      g_nGroup  // Number of groups
    );
  }

  // Main analysis loop: process each marker sequentially
  // Note: OpenMP parallelization is currently disabled due to Rcpp compatibility issues
  for (int i = 0; i < q; i++) {
    // Progress reporting for large analyses
    if (i % 1000 == 0)
      Rcpp::Rcout << "    Completed " << i << "/" << q << " markers in the chunk." << std::endl;

    // Variables to store marker-specific information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;  // Physical position
    bool flip = false;  // Whether a locus is flipped

    uint64_t gIndex = t_genoIndex.at(i);  // Current marker index

    // Extract genotype vector and marker information
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,                      // Genotype file format ("PLINK", "BGEN")
      gIndex,                          // Marker index in genotype file
      ref,                             // Reference allele (output)
      alt,                             // Alternate allele (output)
      marker,                          // Marker ID (output)
      pd,                              // Physical position (output)
      chr,                             // Chromosome (output)
      altFreq,                         // Alternate allele frequency (output)
      altCounts,                       // Alternate allele count (output)
      missingRate,                     // Missing genotype rate (output)
      imputeInfo,                      // Imputation quality score (output)
      true,                            // Whether to output missing data indices
      indexForMissing,                 // Indices of missing genotypes (output)
      false,                           // Whether to output only non-zero genotypes
      indexForNonZero                  // Indices of non-zero genotypes (output)
    );
    int n = GVec.size();  // Sample size

    // Update group-specific statistics if requested
    if (g_ifOutGroup) {
      arma::vec nSamplesInGroupVec(g_nGroup);
      arma::vec AltCountsInGroupVec(g_nGroup);
      arma::vec AltFreqInGroupVec(g_nGroup);

      updateGroupInfo(
        GVec,                     // Genotype vector for current marker
        indexForMissing,          // Indices of individuals with missing genotypes
        nSamplesInGroupVec,       // Output: number of non-missing samples per group
        AltCountsInGroupVec,      // Output: alternate allele counts per group
        AltFreqInGroupVec         // Output: alternate allele frequencies per group
      );

      // Store group statistics in output matrices
      nSamplesInGroup.row(i) = nSamplesInGroupVec.t();
      AltCountsInGroup.row(i) = AltCountsInGroupVec.t();
      AltFreqInGroup.row(i) = AltFreqInGroupVec.t();
    }

    // Format marker information string (CHR:POS:REF:ALT)
    std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;

    // Store basic marker information in output vectors
    markerVec.at(i) = marker;           // Marker ID
    infoVec.at(i) = info;               // Genomic location and alleles
    altFreqVec.at(i) = altFreq;         // Alternate allele frequency
    altCountsVec.at(i) = altCounts;     // Alternate allele count
    missingRateVec.at(i) = missingRate; // Proportion of missing genotypes

    // Calculate MAF and MAC for QC
    double MAF = std::min(
      altFreq,      // Alternate allele frequency
      1 - altFreq   // Reference allele frequency
    );  // MAF is always ≤ 0.5
    double MAC = 2 * MAF * n * (1 - missingRate); // Account for diploid genotypes and missing data

    // Quality Control: check if marker passes MAF, MAC, and missing rate thresholds
    if ((missingRate > g_missingRate_cutoff) ||
        (MAF < g_marker_minMAF_cutoff) ||
        (MAC < g_marker_minMAC_cutoff))
      continue;

    // Check UTIL.cpp
    flip = imputeGenoAndFlip(
      GVec,                // Genotype vector
      altFreq,             // Alternate allele frequency
      indexForMissing,     // Indices of missing genotypes
      missingRate,         // Missing genotype rate
      g_impute_method,     // Imputation method ("mean", "minor", "drop")
      t_method             // Statistical method
    );

    double Beta, seBeta, pval, zScore, hwepval;
    if (t_method == "SPAmix") {
      pval = ptr_gSPAmixobj->getMarkerPval(GVec, altFreq);
      arma::vec pvalVecTemp = ptr_gSPAmixobj->getpvalVec();
      arma::vec zScoreVecTemp = ptr_gSPAmixobj->getzScoreVec();
      
      for (int j = 0; j < Npheno; j++) {
        pvalVec.at(i * Npheno + j) = pvalVecTemp.at(j);
        zScoreVec.at(i * Npheno + j) = zScoreVecTemp.at(j);
      }
      
    } else if (t_method == "SAGELD") {
      pval = ptr_gSAGELDobj->getMarkerPval(GVec, altFreq, hwepval);
      arma::vec pvalVecTemp = ptr_gSAGELDobj->getpvalVec();
      arma::vec zScoreVecTemp = ptr_gSAGELDobj->getzScoreVec();
      arma::vec BetaVecTemp = ptr_gSAGELDobj->getBetaVec();
      arma::vec seBetaVecTemp = ptr_gSAGELDobj->getseBetaVec();

      for (int j = 0; j < 2; j++) {
        pvalVec.at(2 * i + j) = pvalVecTemp.at(j);
        zScoreVec.at(2 * i + j) = zScoreVecTemp.at(j);
        BetaVec.at(2 * i + j) = BetaVecTemp.at(j) * (1 - 2 * flip);
        seBetaVec.at(2 * i + j) = seBetaVecTemp.at(j);
      }    

    } else if (t_method == "WtCoxG") {
      arma::vec pvalVecTemp = ptr_gWtCoxGobj->getpvalVec(GVec, t_genoIndex[i]);
      pvalVec[2 * i]     = pvalVecTemp[0];
      pvalVec[2 * i + 1] = pvalVecTemp[1];

    } else {
      if (t_method == "POLMM") {
        ptr_gPOLMMobj->getMarkerPval(GVec, Beta, seBeta, pval, altFreq, zScore);

      } else if (t_method == "SPACox") {
        pval = ptr_gSPACoxobj->getMarkerPval(GVec, altFreq, zScore);

      } else if (t_method == "SPAGRM") {
        pval = ptr_gSPAGRMobj->getMarkerPval(GVec, altFreq, zScore, hwepval);

      }  else {
        Rcpp::stop("Internal errer in mainMarkerInCPP: Method " + t_method);
      }
      
      pvalVec.at(i) = pval;
      zScoreVec.at(i) = zScore;
      BetaVec.at(i) = Beta * (1 - 2 * flip); // Beta if flip = false, -1*Beta is flip = true
      seBetaVec.at(i) = seBeta;
      hwepvalVec.at(i) = hwepval;
    }
  } // End of one marker loop

  Rcpp::List OutList = Rcpp::List::create(
    Rcpp::Named("markerVec") = markerVec,
    Rcpp::Named("infoVec") = infoVec,
    Rcpp::Named("altFreqVec") = altFreqVec,
    Rcpp::Named("altCountsVec") = altCountsVec,
    Rcpp::Named("missingRateVec") = missingRateVec,
    Rcpp::Named("pvalVec") = pvalVec,
    Rcpp::Named("beta") = BetaVec,
    Rcpp::Named("seBeta") = seBetaVec,
    Rcpp::Named("zScore") = zScoreVec,
    Rcpp::Named("nSamplesInGroup") = nSamplesInGroup,
    Rcpp::Named("AltCountsInGroup") = AltCountsInGroup,
    Rcpp::Named("AltFreqInGroup") = AltFreqInGroup,
    Rcpp::Named("hwepvalVec") = hwepvalVec
  );

  return OutList;
}

//==============================================================================
// SECTION 6: REGION-BASED ANALYSIS FUNCTIONS
// Functions for gene-based and rare variant association testing
//==============================================================================

// Analyzes ultra-rare variants by collapsing them into a single burden score.
// Used for rare variant association analysis where individual variants have
// very low minor allele counts.
// [[Rcpp::export]]
Rcpp::List mainRegionURVInCPP(
  std::string t_method,               // Statistical method for analysis
  std::string t_genoType,             // Genotype file format
  std::vector<uint64_t> t_genoIndex,  // Vector of marker indices in the region
  unsigned int t_n                    // Sample size
) {
  unsigned int q = t_genoIndex.size(); // Number of URV markers after QC
  double Stat, Beta, seBeta, pval0, pval1;
  arma::vec P1Vec(t_n), P2Vec(t_n);

  // Initialize aggregated URV genotype vector (max of individual variants)
  arma::vec GVecURV(t_n, arma::fill::zeros);

  // Process each ultra-rare variant
  for (unsigned int i = 0; i < q; i++) {
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;

    uint64_t gIndex = t_genoIndex.at(i);

    // Extract genotype data for current marker
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,                      // Genotype file format
      gIndex,                          // Marker index in genotype file
      ref,                             // Reference allele (output)
      alt,                             // Alternate allele (output)
      marker,                          // Marker ID (output)
      pd,                              // Physical position (output)
      chr,                             // Chromosome (output)
      altFreq,                         // Alternate allele frequency (output)
      altCounts,                       // Alternate allele count (output)
      missingRate,                     // Missing genotype rate (output)
      imputeInfo,                      // Imputation quality score (output)
      true,                            // Whether to output missing data indices
      indexForMissing,                 // Indices of missing genotypes (output)
      false,                           // Whether to output only non-zero genotypes
      indexForNonZero                  // Indices of non-zero genotypes (output)
    );

    // Impute missing genotypes and flip if necessary
    imputeGenoAndFlip(
      GVec,                // Genotype vector
      altFreq,             // Alternate allele frequency
      indexForMissing,     // Indices of missing genotypes
      missingRate,         // Missing genotype rate
      g_impute_method      // Imputation method ("mean", "minor", "drop")
    );

    // Aggregate variants using max operator (burden test approach)
    if (altFreq < 0.5) {
      GVecURV = arma::max(
        GVecURV, // Current aggregated vector
        GVec     // Current variant genotype vector
      );        // Use alternate allele
    } else {
      GVecURV = arma::max(
        GVecURV,  // Current aggregated vector
        2 - GVec  // Flipped genotype vector
      );    // Use reference allele (flip)
    }
  }

  // Perform region-level association test on aggregated genotype
  Unified_getRegionPVec(
    t_method,                 // Statistical method
    GVecURV,                  // Aggregated URV genotype vector
    Stat,                     // Score statistic (output)
    Beta,                     // Effect size estimate (output)
    seBeta,                   // Standard error of effect size (output)
    pval0,                    // Normal approximation p-value (output)
    pval1,                    // SPA-corrected p-value (output)
    P1Vec,                    // First component vector for variance calculation (output)
    P2Vec                     // Second component vector for variance calculation (output)
  );

  // Return results
  Rcpp::List OutList = Rcpp::List::create(
    Rcpp::Named("Stat") = Stat,
    Rcpp::Named("pval1") = pval1
  );

  return OutList;
}

// Performs comprehensive gene-based analysis including:
// - Individual variant association tests
// - Burden tests with different MAF thresholds
// - Ultra-rare variant collapsing
// - Annotation-stratified analysis
// - Variance component tests
// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(
  std::string t_method,                      // Statistical method
  std::string t_genoType,                    // Genotype file format
  std::vector<uint64_t> t_genoIndex,         // Vector of marker indices in region
  std::vector<double> t_weightVec,           // Variant weights for burden testing
  std::string t_outputFile,                  // Base filename for temporary output files
  std::vector<unsigned int> t_labelVec,      // Individual labels for stratified analysis
  unsigned int t_nLabel,                     // Number of distinct labels
  arma::mat t_annoMat,                       // Annotation matrix (markers × annotations)
  std::vector<std::string> t_annoVec         // Annotation names
) {
  unsigned int n = t_labelVec.size();                     // Sample size
  unsigned int q = t_genoIndex.size();                    // Number of markers in region
  unsigned int nAnno = t_annoMat.n_cols;                  // Number of annotations

  // Initialize output vectors and matrices
  arma::uvec indicatorVec(
    q + nAnno,        // Total number of markers and annotations
    arma::fill::zeros // Initialize with zeros
  );                  // QC status: 0=fail, 1=pass, 2=URV, 3=aggregated
  Rcpp::StringVector markerVec(q + nAnno);                 // Marker IDs
  Rcpp::StringVector infoVec(q + nAnno);                   // Marker information
  arma::vec altFreqVec(q + nAnno);                         // Alternate allele frequencies
  arma::vec MACVec(q + nAnno);                             // Minor allele counts
  arma::vec MAFVec(q + nAnno);                             // Minor allele frequencies
  arma::vec missingRateVec(q + nAnno);                     // Missing rates

  // Statistical results vectors
  std::vector<double> altBetaVec(q + nAnno);               // Effect sizes
  std::vector<double> seBetaVec(q + nAnno);                // Standard errors
  std::vector<double> pval0Vec(q + nAnno);                 // Normal approximation p-values
  std::vector<double> pval1Vec(q + nAnno);                 // SPA-corrected p-values
  std::vector<double> StatVec(q + nAnno);                  // Score statistics

  // Label-stratified analysis matrices
  arma::mat MACLabelMat(
    q + nAnno,  // Total markers and annotations
    t_nLabel    // Number of label categories
  );            // MAC by label
  arma::mat MAFLabelMat(
    q + nAnno,  // Total markers and annotations
    t_nLabel    // Number of label categories
  );            // MAF by label

  // Memory management for large regions
  unsigned int m1 = g_region_maxMarkers_cutoff;           // Markers per chunk
  arma::mat P1Mat(
    m1,  // Maximum markers per chunk
    n    // Number of samples
  );     // Component matrix 1
  arma::mat P2Mat(
    n,   // Number of samples
    m1   // Maximum markers per chunk
  );     // Component matrix 2
  std::vector<unsigned int> mPassCVVec;                   // Markers passing QC per chunk

  // Analysis state variables
  unsigned int nchunks = 0;       // Total number of chunks
  unsigned int ichunk = 0;        // Current chunk index
  unsigned int i1InChunk = 0;     // Non-URV markers in current chunk
  unsigned int i1 = 0;            // Total non-URV markers processed
  unsigned int i2 = 0;            // Total URV markers processed

  // Aggregation matrices for different analysis types
  arma::mat GMatURV(
    n,                // Number of samples
    nAnno,            // Number of annotations
    arma::fill::zeros // Initialize with zeros
  );                  // Ultra-rare variant matrix
  unsigned int n_max_maf = g_region_max_maf_vec.size();    // Number of MAF thresholds
  arma::mat GMatBurden(
    n,                     // Number of samples
    nAnno * n_max_maf,     // Number of burden test combinations
    arma::fill::zeros      // Initialize with zeros
  );                       // Weighted burden matrix
  arma::mat pvalBurden(
    nAnno * n_max_maf,     // Number of burden test combinations
    2,                     // Two p-values per test (normal, SPA)
    arma::fill::zeros      // Initialize with zeros
  );                       // Burden test p-values
  arma::mat GMatBurdenNoWeight(
    n,                     // Number of samples
    nAnno * n_max_maf,     // Number of burden test combinations
    arma::fill::zeros      // Initialize with zeros
  );                       // Unweighted burden matrix
  arma::mat infoBurdenNoWeight(
    nAnno * n_max_maf,     // Number of burden test combinations
    7,                     // Number of information columns
    arma::fill::zeros      // Initialize with zeros
  );                       // Burden test information

  // Beta distribution for variant weighting
  boost::math::beta_distribution<> beta_dist(
    g_region_weight_beta[0], // Alpha parameter for beta distribution
    g_region_weight_beta[1]  // Beta parameter for beta distribution
  );

  // MAIN ANALYSIS LOOP: Process each marker in the region
  for (unsigned int i = 0; i < q; i++) {
    double weight = t_weightVec.at(i);                 // Variant-specific weight

    // Extract marker information and genotype data
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    uint64_t gIndex = t_genoIndex.at(i);

    // Timing analysis for performance monitoring
    arma::vec test11 = getTime();

    arma::vec GVec = Unified_getOneMarker(
      t_genoType,                      // Genotype file format
      gIndex,                          // Marker index in genotype file
      ref,                             // Reference allele (output)
      alt,                             // Alternate allele (output)
      marker,                          // Marker ID (output)
      pd,                              // Physical position (output)
      chr,                             // Chromosome (output)
      altFreq,                         // Alternate allele frequency (output)
      altCounts,                       // Alternate allele count (output)
      missingRate,                     // Missing genotype rate (output)
      imputeInfo,                      // Imputation quality score (output)
      true,                            // Whether to output missing data indices
      indexForMissing,                 // Indices of missing genotypes (output)
      false,                           // Whether to output only non-zero genotypes
      indexForNonZero                  // Indices of non-zero genotypes (output)
    );

    arma::vec test12 = getTime();
    g_compTime1 += test12 - test11;

    // Format marker information
    std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;

    // Impute missing genotypes and determine allele flipping
    flip = imputeGenoAndFlip(
      GVec,                // Genotype vector
      altFreq,             // Alternate allele frequency
      indexForMissing,     // Indices of missing genotypes
      missingRate,         // Missing genotype rate
      g_impute_method      // Imputation method ("mean", "minor", "drop")
    );

    // Calculate quality control metrics
    double MAF = std::min(
      altFreq,      // Alternate allele frequency
      1 - altFreq   // Reference allele frequency
    );
    double MAC = MAF * 2 * n * (1 - missingRate);

    // Store basic marker information
    markerVec.at(i) = marker;
    infoVec.at(i) = info;
    altFreqVec.at(i) = altFreq;
    missingRateVec.at(i) = missingRate;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;

    // Quality Control: Skip markers that don't meet criteria
    if ((missingRate > g_missingRate_cutoff) || (MAF > g_region_maxMAF_cutoff) || MAF == 0) {
      continue; // Marker fails QC
    }

    if (MAC > g_region_minMAC_cutoff) {
      // NON-ULTRA-RARE VARIANT: Perform individual marker analysis
      indicatorVec.at(i) = 1;

      if (i1InChunk == 0) {
        Rcpp::Rcout << "    Start analyzing chunk " << ichunk << " ..." << std::endl;
      }

      // Statistical analysis timing
      arma::vec test21 = getTime();

      // Individual marker association testing
      double Stat, Beta, seBeta, pval0, pval1;
      arma::vec P1Vec(n), P2Vec(n);

      Unified_getRegionPVec(
        t_method,             // Statistical method
        GVec,                 // Genotype vector for current marker
        Stat,                 // Score statistic (output)
        Beta,                 // Effect size estimate (output)
        seBeta,               // Standard error of effect size (output)
        pval0,                // Normal approximation p-value (output)
        pval1,                // SPA-corrected p-value (output)
        P1Vec,                // First component vector for variance calculation (output)
        P2Vec                 // Second component vector for variance calculation (output)
      );

      arma::vec test22 = getTime();
      g_compTime2 += test22 - test21;

      // Store statistical results
      StatVec.at(i) = Stat;
      altBetaVec.at(i) = Beta * (1 - 2 * flip);  // Adjust for allele flipping
      seBetaVec.at(i) = seBeta;
      pval0Vec.at(i) = pval0;
      pval1Vec.at(i) = pval1;

      // Store variance components for region-level testing
      P1Mat.row(i1InChunk) = P1Vec.t();
      P2Mat.col(i1InChunk) = P2Vec;

      // Label-stratified analysis
      if (t_nLabel != 1)
        getLabelInfo(
          GVec,           // Genotype vector after imputation
          t_labelVec,     // Label assignment for each individual
          i,              // Row index in output matrices
          MACLabelMat,    // Output matrix for Minor Allele Counts by label
          MAFLabelMat     // Output matrix for Minor Allele Frequencies by label
        );

      i1 += 1;
      i1InChunk += 1;

      // Burden test aggregation with beta-function weighting
      double w0 = boost::math::pdf(
        beta_dist, // Beta distribution object
        MAF        // Minor allele frequency value
      );  // Beta-distribution weight

      for (unsigned int j = 0; j < n; j++) {
        if (GVec.at(j) != 0) {
          for (unsigned int iAnno = 0; iAnno < nAnno; iAnno++) {
            if (t_annoMat(i, iAnno) == 1)  // Marker in this annotation
            {
              for (unsigned int i_max_maf = 0; i_max_maf < n_max_maf; i_max_maf++) {
                double max_maf = g_region_max_maf_vec.at(i_max_maf);
                if (MAF < max_maf) {
                  // Add to weighted and unweighted burden scores
                  GMatBurden(j, iAnno * n_max_maf + i_max_maf) += w0 * GVec.at(j);
                  GMatBurdenNoWeight(j, iAnno * n_max_maf + i_max_maf) += GVec.at(j);
                }
              }
            }
          }
        }
      }
    } else {
      // ULTRA-RARE VARIANT: Aggregate for burden testing
      indicatorVec.at(i) = 2;

      for (unsigned int j = 0; j < n; j++) {
        if (GVec.at(j) != 0) {
          for (unsigned iAnno = 0; iAnno < nAnno; iAnno++) {
            if (t_annoMat(i, iAnno) == 1) {
              // Use max aggregation for ultra-rare variants
              GMatURV(j, iAnno) = std::max(
                GMatURV(j, iAnno),    // Current aggregated value
                weight * GVec.at(j)   // Weighted current variant value
              );

              // Also add to burden scores
              for (unsigned int i_max_maf = 0; i_max_maf < n_max_maf; i_max_maf++) {
                GMatBurdenNoWeight(j, iAnno * n_max_maf + i_max_maf) += GVec.at(j);
              }
            }
          }
        }
      }
      i2 += 1;
    }

    // Chunk management: Save intermediate results if chunk is full
    if (i1InChunk == m1) {
      Rcpp::Rcout << "    In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and "
                << i1 << " markers are not ultra-rare." << std::endl;

      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");

      mPassCVVec.push_back(m1);
      ichunk += 1;
      i1InChunk = 0;
    }

    Rcpp::checkUserInterrupt();  // Allow R to interrupt long-running analysis
  }

  // Print timing information
  printTimeDiff(
    g_compTime1,                // Timing vector for function
    "Unified_getOneMarker"      // Function name description
  );
  printTimeDiff(
    g_compTime2,                // Timing vector for function
    "Unified_getRegionPVec"     // Function name description
  );

  // ULTRA-RARE VARIANT ANALYSIS: Process aggregated URVs by annotation
  for (unsigned int iAnno = 0; iAnno < nAnno; iAnno++) {
    arma::vec GVecURV = GMatURV.col(iAnno);

    // Burden testing with beta-function weighting for URVs
    double MAFURV = mean(GVecURV) / 2;
    double w0URV = boost::math::pdf(
      beta_dist, // Beta distribution object
      MAFURV     // Ultra-rare variant frequency value
    );

    for (unsigned int i_max_maf = 0; i_max_maf < n_max_maf; i_max_maf++) {
      unsigned int i_pos = iAnno * n_max_maf + i_max_maf;

      // Add weighted URV contribution to burden scores
      GMatBurden.col(i_pos) += w0URV * GVecURV;

      // Test weighted burden score
      double Stat, Beta, seBeta, pval0, pval1;
      arma::vec P1Vec(n), P2Vec(n);

      Unified_getRegionPVec(
        t_method,                     // Statistical method
        GMatBurden.col(i_pos),        // Weighted burden vector for this annotation/MAF combination
        Stat,                         // Score statistic (output)
        Beta,                         // Effect size estimate (output)
        seBeta,                       // Standard error of effect size (output)
        pval0,                        // Normal approximation p-value (output)
        pval1,                        // SPA-corrected p-value (output)
        P1Vec,                        // First component vector for variance calculation (output)
        P2Vec                         // Second component vector for variance calculation (output)
      );
      pvalBurden.at(i_pos, 0) = pval0;
      pvalBurden.at(i_pos, 1) = pval1;

      // Test unweighted burden score and store detailed information
      Unified_getRegionPVec(
        t_method,                     // Statistical method
        GMatBurdenNoWeight.col(i_pos), // Unweighted burden vector for this annotation/MAF combination
        Stat,                         // Score statistic (output)
        Beta,                         // Effect size estimate (output)
        seBeta,                       // Standard error of effect size (output)
        pval0,                        // Normal approximation p-value (output)
        pval1,                        // SPA-corrected p-value (output)
        P1Vec,                        // First component vector for variance calculation (output)
        P2Vec                         // Second component vector for variance calculation (output)
      );
      infoBurdenNoWeight.at(i_pos, 0) = iAnno;                              // Annotation index
      infoBurdenNoWeight.at(i_pos, 1) = i_max_maf;                          // MAF threshold index
      infoBurdenNoWeight.at(i_pos, 2) = sum(GMatBurdenNoWeight.col(i_pos)); // Total burden score
      infoBurdenNoWeight.at(i_pos, 3) = Stat;                               // Score statistic
      infoBurdenNoWeight.at(i_pos, 4) = Beta;                               // Effect size
      infoBurdenNoWeight.at(i_pos, 5) = seBeta;                             // Standard error
      infoBurdenNoWeight.at(i_pos, 6) = pval1;                              // SPA p-value
    }

    // Store URV aggregation results
    indicatorVec.at(q + iAnno) = 3;                           // Aggregated URV indicator
    markerVec.at(q + iAnno) = t_annoVec.at(iAnno);            // Annotation name
    infoVec.at(q + iAnno) = "Ultra-Rare Variants";            // Description
    altFreqVec.at(q + iAnno) = MAFVec.at(q + iAnno) = MAFURV; // Aggregated frequency
    MACVec.at(q + iAnno) = sum(GVecURV);                      // Total allele count
    missingRateVec.at(q + iAnno) = 0;                         // No missing data after aggregation

    // Label-stratified analysis for URV
    if (t_nLabel != 1)
      getLabelInfo(
        GVecURV,        // Genotype vector after imputation
        t_labelVec,     // Label assignment for each individual
        q + iAnno,      // Row index in output matrices
        MACLabelMat,    // Output matrix for Minor Allele Counts by label
        MAFLabelMat     // Output matrix for Minor Allele Frequencies by label
      );

    // Statistical testing for aggregated URV
    double Stat, Beta, seBeta, pval0, pval1;
    arma::vec P1Vec(n), P2Vec(n);

    Unified_getRegionPVec(
      t_method,  // Statistical method
      GVecURV,   // Aggregated ultra-rare variant genotype vector
      Stat,      // Score statistic (output)
      Beta,      // Effect size estimate (output)
      seBeta,    // Standard error of effect size (output)
      pval0,     // Normal approximation p-value (output)
      pval1,     // SPA-corrected p-value (output)
      P1Vec,     // First component vector for variance calculation (output)
      P2Vec      // Second component vector for variance calculation (output)
    );

    // Store URV statistical results
    StatVec.at(q + iAnno) = Stat;
    altBetaVec.at(q + iAnno) = Beta;
    seBetaVec.at(q + iAnno) = seBeta;
    pval0Vec.at(q + iAnno) = pval0;
    pval1Vec.at(q + iAnno) = pval1;

    // Add URV results to variance component matrices
    if (i1InChunk >= m1) {
      P1Mat.resize(
        i1InChunk + 1,  // New number of rows
        n               // Number of columns (samples)
      );
      P2Mat.resize(
        n,              // Number of rows (samples)
        i1InChunk + 1   // New number of columns
      );
    }

    P1Mat.row(i1InChunk) = P1Vec.t();
    P2Mat.col(i1InChunk) = P2Vec;
    i1 += 1;
    i1InChunk += 1;
  }

  // Validation checks
  if (i2 == 0)
    Rcpp::Rcout << "    i2 == 0." << std::endl;

  if ((i1 == 0) & (i2 == 0))
    Rcpp::Rcout << "    Cannot find any valid rare variants. This region will be skipped." << std::endl;

  // Finalize chunk management
  mPassCVVec.push_back(i1InChunk);
  nchunks = ichunk + 1;

  // Save final chunk if needed
  if (i1InChunk != 0) {
    P1Mat = P1Mat.rows(
      0,              // Starting row
      i1InChunk - 1   // Ending row
    );
    P2Mat = P2Mat.cols(
      0,              // Starting column
      i1InChunk - 1   // Ending column
    );

    if (nchunks != 1) {
      Rcpp::Rcout << "    In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and "
                << i1 << " markers are not ultra-rare." << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
  }

  // VARIANCE MATRIX COMPUTATION for region-level testing
  arma::mat VarMat(i1, i1);

  if (nchunks == 1) {
    // All data fits in memory
    VarMat = P1Mat * P2Mat;
  } else {
    // Compute variance matrix from chunks stored on disk
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;

    for (unsigned int index1 = 0; index1 < nchunks; index1++) {
      last_row = first_row + mPassCVVec.at(index1) - 1;
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      P1Mat.load(P1MatFile);

      if (P1Mat.n_cols == 0)
        continue;

      // Compute off-diagonal blocks
      for (unsigned int index2 = 0; index2 < index1; index2++) {
        Rcpp::Rcout << "    Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", "
                  << index2 << "/" << nchunks - 1 << ") ..." << std::endl;

        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");

        if (P2Mat.n_cols == 0)
          continue;

        arma::mat offVarMat = P1Mat * P2Mat;
        last_col = first_col + mPassCVVec.at(index2) - 1;

        VarMat.submat(
          first_row,  // Starting row
          first_col,  // Starting column
          last_row,   // Ending row
          last_col    // Ending column
        ) = offVarMat;
        VarMat.submat(
          first_col,  // Starting row (transposed)
          first_row,  // Starting column (transposed)
          last_col,   // Ending row (transposed)
          last_row    // Ending column (transposed)
        ) = offVarMat.t(); // Symmetric

        first_col = last_col + 1;
      }

      // Compute diagonal block
      last_col = first_col + mPassCVVec.at(index1) - 1;
      Rcpp::Rcout << "    Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", "
                << index1 << "/" << nchunks - 1 << ") ..." << std::endl;

      P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
      arma::mat diagVarMat = P1Mat * P2Mat;
      VarMat.submat(
        first_row,  // Starting row
        first_col,  // Starting column
        last_row,   // Ending row
        last_col    // Ending column
      ) = diagVarMat;

      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
    }

    // Clean up temporary files
    for (unsigned int index1 = 0; index1 < nchunks; index1++) {
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      std::string P2MatFile = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin";
      const char *File1 = P1MatFile.c_str();
      const char *File2 = P2MatFile.c_str();
      std::remove(File1);
      std::remove(File2);
    }
  }

  // Rcpp::Rcout << "    m1:\t" << m1 << std::endl;

  // PREPARE OUTPUT: Comprehensive results list
  Rcpp::List OutList = Rcpp::List::create(
    Rcpp::Named("markerVec") = markerVec,                        // Marker IDs
    Rcpp::Named("infoVec") = infoVec,                            // Marker information
    Rcpp::Named("missingRateVec") = missingRateVec,              // Missing rates
    Rcpp::Named("altFreqVec") = altFreqVec,                      // Allele frequencies
    Rcpp::Named("MACVec") = MACVec,                              // Minor allele counts
    Rcpp::Named("MAFVec") = MAFVec,                              // Minor allele frequencies
    Rcpp::Named("MACLabelMat") = MACLabelMat,                    // MAC by label
    Rcpp::Named("MAFLabelMat") = MAFLabelMat,                    // MAF by label
    Rcpp::Named("StatVec") = StatVec,                            // Score statistics
    Rcpp::Named("altBetaVec") = altBetaVec,                      // Effect sizes
    Rcpp::Named("seBetaVec") = seBetaVec,                        // Standard errors
    Rcpp::Named("pval0Vec") = pval0Vec,                          // Normal p-values
    Rcpp::Named("pval1Vec") = pval1Vec,                          // SPA p-values
    Rcpp::Named("indicatorVec") = indicatorVec,                  // QC indicators
    Rcpp::Named("VarMat") = VarMat,                              // Variance matrix
    Rcpp::Named("pvalBurden") = pvalBurden,                      // Burden p-values
    Rcpp::Named("infoBurdenNoWeight") = infoBurdenNoWeight       // Unweighted burden info
  );

  return OutList;
}


//==============================================================================
// SECTION 8: GENOTYPE EXTRACTION AND MATRIX FUNCTIONS
// Functions for extracting and formatting genotype data
//==============================================================================

// Efficiently extracts basic genotype statistics without loading full genotype data.
// Useful for quality control and marker filtering before analysis.
// [[Rcpp::export]]
arma::mat getGenoInfoInCPP(
  std::string t_genoType,       // Genotype file format
  Rcpp::DataFrame t_markerInfo  // DataFrame containing marker information and indices
) {
  int q = t_markerInfo.nrow();  // number of markers requested
  arma::mat genoInfoMat(q, 2);

  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];

  // Variables for storing marker information from Unified_getOneMarker
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::string ref, alt, marker, chr;
  std::vector<uint32_t> indexForMissing, indexForNonZero;

  for (int i = 0; i < q; i++) {
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,           // Genotype file format ("PLINK", "BGEN")
      gIndex,               // Marker index (different meanings for different formats)
      ref,                  // Reference allele (output)
      alt,                  // Alternate allele (output, should be minor allele for efficiency)
      marker,               // Marker ID extracted from genotype file (output)
      pd,                   // Base position (output)
      chr,                  // Chromosome (output)
      altFreq,              // Frequency of alternate allele (output)
      altCounts,            // Counts of alternate allele (output)
      missingRate,          // Missing rate (output)
      imputeInfo,           // Imputation information score, R2 (all 1 for PLINK) (output)
      true,                 // If true, output index of missing genotype data
      indexForMissing,      // Index of missing genotype data (output)
      false,                // If true, only output non-zero genotypes (NOTE: inefficient if ALT isn't minor)
      indexForNonZero       // Index of non-zero genotypes (output, only valid if above is true)
    );

    if ((i + 1) % 1000 == 0)
      Rcpp::Rcout << "    Completed " << (i + 1) << "/" << q << " genetic variants." << std::endl;

    genoInfoMat.at(i, 0) = altFreq;
    genoInfoMat.at(i, 1) = missingRate;
  }

  return genoInfoMat;
}

// [[Rcpp::export]]
arma::mat getGenoInCPP(
  std::string t_genoType,
  Rcpp::DataFrame t_markerInfo,
  int n,
  std::string t_imputeMethod  // "none", "mean", "bestguess"
) {
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::mat GMat(n, q);

  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;

  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,          // "PLINK", "BGEN"
      gIndex,              // different meanings for different genoType
      ref,                 // REF allele
      alt,                 // ALT allele (should probably be minor allele, otherwise, computation time will increase)
      marker,              // marker ID extracted from genotype file
      pd,                  // base position
      chr,                 // chromosome
      altFreq,             // frequency of ALT allele
      altCounts,           // counts of ALT allele
      missingRate,         // missing rate
      imputeInfo,          // imputation information score, i.e., R2 (all 1 for PLINK)
      true,                // if true, output index of missing genotype data
      indexForMissing,     // index of missing genotype data
      false,               // if true, only output a vector of non-zero genotype
      indexForNonZero      // the index of non-zero genotype. Only valid if t_isOnlyOutputNonZero == true
    );
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.cpp
    GMat.col(i) = GVec;
  }

  return GMat;
}

// This function will replace the above function (2022-01-28)
// [[Rcpp::export]]
arma::mat getGenoInCPP_fixedNumber(
  std::string t_genoType,        // Genotype file format
  Rcpp::DataFrame t_markerInfo,  // DataFrame containing marker information
  int n,                         // Sample size
  std::string t_imputeMethod,    // "none", "mean", "bestguess"
  int m,                         // Number of selected markers
  double missingRateCutoff,      // Maximum allowed missing rate
  double minMAFCutoff            // Minimum MAF threshold
) {
  int q = t_markerInfo.nrow(); // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];

  arma::mat GMat(n, m);
  int index = 0;

  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;

  for (int i = 0; i < q; i++) {
    uint64_t gIndex = gIndexVec.at(i);
    // Rcpp::Rcout << "    gIndex:\t" << gIndex << std::endl;
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,             // Genotype file format ("PLINK", "BGEN")
      gIndex,                 // Marker index (different meanings for different formats)
      ref,                    // Reference allele (output)
      alt,                    // Alternate allele (output, should be minor allele for efficiency)
      marker,                 // Marker ID extracted from genotype file (output)
      pd,                     // Base position (output)
      chr,                    // Chromosome (output)
      altFreq,                // Frequency of alternate allele (output)
      altCounts,              // Counts of alternate allele (output)
      missingRate,            // Missing rate (output)
      imputeInfo,             // Imputation information score, R2 (all 1 for PLINK) (output)
      true,                   // If true, output index of missing genotype data
      indexForMissing,        // Index of missing genotype data (output)
      false,                  // If true, only output non-zero genotypes (NOTE: inefficient if ALT isn't minor)
      indexForNonZero         // Index of non-zero genotypes (output, only valid if t_isOnlyOutputNonZero == true)
    );

    if ((altFreq < minMAFCutoff) | (altFreq > 1 - minMAFCutoff) | (missingRate > missingRateCutoff))
      continue;

    imputeGeno(
      GVec,               // Genotype vector
      altFreq,            // Alternate allele frequency
      indexForMissing,    // Indices of missing genotypes
      t_imputeMethod      // Imputation method
    ); // check UTIL.cpp
    GMat.col(index) = GVec;
    index++;
    if (index == m)
      break;
  }

  if (index < m)
    Rcpp::stop("No enough variants are for variance ratio estimation.");

  return GMat;
}

// [[Rcpp::export]]
arma::sp_mat getSpGenoInCPP(
  std::string t_genoType,       // Genotype file format
  Rcpp::DataFrame t_markerInfo, // DataFrame containing marker information
  int n,                        // Sample size
  std::string t_imputeMethod    // Imputation method
) {
  int q = t_markerInfo.nrow(); // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::sp_mat GMat(n, q); // change #1 compared to getGenoInCPP()

  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;

  for (int i = 0; i < q; i++) {
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(
      t_genoType,             // Genotype file format ("PLINK", "BGEN")
      gIndex,                 // Marker index (different meanings for different formats)
      ref,                    // Reference allele (output)
      alt,                    // Alternate allele (output, should be minor allele for efficiency)
      marker,                 // Marker ID extracted from genotype file (output)
      pd,                     // Base position (output)
      chr,                    // Chromosome (output)
      altFreq,                // Frequency of alternate allele (output)
      altCounts,              // Counts of alternate allele (output)
      missingRate,            // Missing rate (output)
      imputeInfo,             // Imputation information score, R2 (all 1 for PLINK) (output)
      true,                   // If true, output index of missing genotype data
      indexForMissing,        // Index of missing genotype data (output)
      false,                  // If true, only output non-zero genotypes (NOTE: inefficient if ALT isn't minor)
      indexForNonZero         // Index of non-zero genotypes (output, only valid if above is true)
    );

    imputeGeno(
      GVec,               // Genotype vector
      altFreq,            // Alternate allele frequency
      indexForMissing,    // Indices of missing genotypes
      t_imputeMethod      // Imputation method
    );
    GMat.col(i) = arma::sp_mat(GVec);            // change #2 compared to getGenoInCPP()
  }

  return GMat;
}


//==============================================================================
// SECTION 9: OBJECT INITIALIZATION AND CLEANUP FUNCTIONS
// Functions for setting up and managing various analysis objects
//==============================================================================

// Creates and configures a PLINK file reader for efficient genotype data access.
// The PLINK format is widely used for large-scale genetic data storage.
// [[Rcpp::export]]
void setPLINKobjInCPP(
  std::string t_bimFile,                    // Path to PLINK .bim file (marker information)
  std::string t_famFile,                    // Path to PLINK .fam file (sample information)
  std::string t_bedFile,                    // Path to PLINK .bed file (binary genotype data)
  std::vector<std::string> t_SampleInModel, // Vector of sample IDs to include in analysis
  std::string t_AlleleOrder                 // Allele ordering convention ("alt-first" or "ref-first")
) {
  if (ptr_gPLINKobj)
    delete ptr_gPLINKobj;

  ptr_gPLINKobj = new PLINK::PlinkClass(
    t_bimFile,                              // Path to PLINK .bim file (marker information)
    t_famFile,                              // Path to PLINK .fam file (sample information)
    t_bedFile,                              // Path to PLINK .bed file (binary genotype data)
    t_SampleInModel,                        // Vector of sample IDs to include in analysis
    t_AlleleOrder                           // Allele ordering convention ("alt-first" or "ref-first")
  );

  int n = ptr_gPLINKobj->getN();
  Rcpp::Rcout << "    Number of subjects with genotype: " << n << std::endl;
}

// Creates and configures a BGEN file reader for efficient genotype data access.
// BGEN is a compressed binary format for large genetic datasets with probabilistic genotypes.
// [[Rcpp::export]]
void setBGENobjInCPP(
  std::string t_bgenFileName,                // Path to BGEN file containing genotype data
  std::string t_bgenFileIndex,               // Path to BGEN index file for fast access
  std::vector<std::string> t_SampleInBgen,   // Sample IDs as they appear in BGEN file
  std::vector<std::string> t_SampleInModel,  // Sample IDs to include in analysis
  bool t_isSparseDosageInBgen,               // Whether BGEN file uses sparse dosage encoding
  bool t_isDropmissingdosagesInBgen,         // Whether to drop missing dosage values
  std::string t_AlleleOrder                  // Allele ordering convention ("alt-first" or "ref-first")
) {
  if (ptr_gBGENobj) {
    // Rcpp::Rcout << "    Deleting `ptr_gBGENobj` ..." << std::endl;
    delete ptr_gBGENobj;
    ptr_gBGENobj = nullptr;
  }

  ptr_gBGENobj = new BGEN::BgenClass(
    t_bgenFileName,                          // Path to BGEN file containing genotype data
    t_bgenFileIndex,                         // Path to BGEN index file for fast access
    t_SampleInBgen,                          // Sample IDs as they appear in BGEN file
    t_SampleInModel,                         // Sample IDs to include in analysis
    t_isSparseDosageInBgen,                  // Whether BGEN file uses sparse dosage encoding
    t_isDropmissingdosagesInBgen,            // Whether to drop missing dosage values
    t_AlleleOrder                            // Allele ordering convention ("alt-first" or "ref-first")
  );
  int n = ptr_gBGENobj->getN();
  Rcpp::Rcout << "    Number of subjects with genotype: " << n << std::endl;
}

// Safely closes and cleans up genotype file readers to prevent memory leaks.
// Should be called when switching between file formats or at end of analysis.
// [[Rcpp::export]]
void closeGenoInputInCPP(
  std::string t_genoType  // Genotype file format to close ("PLINK" or "BGEN")
) {
  if (t_genoType == "PLINK") {
    delete ptr_gPLINKobj;
    ptr_gPLINKobj = nullptr;
  }
  if (t_genoType == "BGEN") {
    delete ptr_gBGENobj;
    ptr_gBGENobj = nullptr;
  }
}

//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

//==============================================================================
// SECTION 10: STATISTICAL METHOD OBJECT INITIALIZATION
// Functions for setting up different statistical analysis methods
//==============================================================================

// Initialize POLMM object for polytomous trait analysis with pre-computed matrices.
// This version uses already computed mean and correlation matrices from prior fitting.
// [[Rcpp::export]]
void setPOLMMobjInCPP(
  arma::mat t_muMat,                         // Mean matrix for logistic mixed model
  arma::mat t_iRMat,                         // Inverse correlation matrix
  arma::mat t_Cova,                          // Covariate matrix for fixed effects
  arma::uvec t_yVec,                         // Response vector (categorical outcomes)
  double t_tau,                              // Variance component parameter
  bool t_printPCGInfo,                       // Whether to print PCG solver information
  double t_tolPCG,                           // Tolerance for PCG solver convergence
  int t_maxiterPCG,                          // Maximum iterations for PCG solver
  double t_varRatio,                         // Variance ratio estimate
  double t_SPA_cutoff,                       // P-value cutoff for SPA correction
  bool t_flagSparseGRM,                      // Whether to use sparse or dense GRM
  arma::uvec t_group,                        // Group assignment for each individual
  bool t_ifOutGroup,                         // Whether to output group-specific statistics
  unsigned int t_nGroup                      // Total number of groups
) {
  // Set POLMM-specific group variables
  g_group = t_group;
  g_ifOutGroup = t_ifOutGroup;
  g_nGroup = t_nGroup;

  if (ptr_gPOLMMobj)
    delete ptr_gPOLMMobj;

  ptr_gPOLMMobj = new POLMM::POLMMClass(
    t_muMat,                                 // Mean matrix
    t_iRMat,                                 // Inverse correlation matrix
    t_Cova,                                  // Covariate matrix
    t_yVec,                                  // Response vector
    g_SparseGRM,                             // Global sparse GRM
    t_tau,                                   // Variance component
    t_printPCGInfo,                          // PCG info printing flag
    t_tolPCG,                                // PCG tolerance
    t_maxiterPCG,                            // PCG max iterations
    t_varRatio,                              // Variance ratio
    t_SPA_cutoff,                            // SPA cutoff threshold
    t_flagSparseGRM                          // Sparse GRM flag
  );
}

// Initialize and fit POLMM null model for polytomous traits with optional LOCO analysis.
// This version handles null model fitting and variance ratio estimation.
// [[Rcpp::export]]
Rcpp::List setPOLMMobjInCPP_NULL(
  bool t_flagSparseGRM,                      // If true, use sparse GRM; otherwise use dense GRM
  arma::mat t_Cova,                          // Covariate matrix for fixed effects
  arma::uvec t_yVec,                         // Response vector (values from 0 to J-1 for J categories)
  arma::vec t_beta,                          // Fixed effect coefficients
  arma::vec t_bVec,                          // Random effect coefficients
  arma::vec t_eps,                           // Residual vector
  double t_tau,                              // Variance component parameter
  Rcpp::List t_SPmatR,                       // Sparse matrix representation (output of makeSPmatR())
  Rcpp::List t_controlList,                  // Control parameters for optimization
  arma::mat GenoMat                          // Genotype matrix for variance ratio estimation
) {
  if (ptr_gPOLMMobj)
    delete ptr_gPOLMMobj;

  ptr_gPOLMMobj = new POLMM::POLMMClass(
    t_flagSparseGRM,                         // If true, use sparse GRM; otherwise use dense GRM
    ptr_gDenseGRMobj,                        // Dense GRM object pointer
    ptr_gPLINKobj,                           // PLINK file reader object pointer
    ptr_gBGENobj,                            // BGEN file reader object pointer
    t_Cova,                                  // Covariate matrix for fixed effects
    t_yVec,                                  // Response vector (values from 0 to J-1)
    t_beta,                                  // Fixed effect coefficients
    t_bVec,                                  // Random effect coefficients
    t_eps,                                   // Residual vector
    t_tau,                                   // Variance component parameter
    g_SparseGRM,                             // Global sparse GRM (from setSparseGRMInCPP())
    t_controlList                            // Control parameters for optimization
  );

  ptr_gPOLMMobj->fitPOLMM();
  ptr_gPOLMMobj->estVarRatio(GenoMat);

  Rcpp::List outList = ptr_gPOLMMobj->getPOLMM();
  return outList;
}

// Initialize SPAGRM object for saddle point approximation with genetic relationship matrix.
// Handles outlier-robust analysis with pre-computed residuals and GRM components.
// [[Rcpp::export]]
void setSPAGRMobjInCPP(
  arma::vec t_resid,                         // Residual vector from null model
  arma::vec t_resid_unrelated_outliers,      // Residuals for unrelated outlier subjects
  double t_sum_R_nonOutlier,                 // Sum of residuals for non-outlier subjects
  double t_R_GRM_R_nonOutlier,               // Quadratic form R'*GRM*R for non-outliers
  double t_R_GRM_R_TwoSubjOutlier,           // Quadratic form for two-subject outlier pairs
  double t_R_GRM_R,                          // Full quadratic form R'*GRM*R
  arma::vec t_MAF_interval,                  // Minor allele frequency intervals for binning
  Rcpp::List t_TwoSubj_list,                 // List of two-subject outlier pairs information
  Rcpp::List t_ThreeSubj_list,               // List of three-subject outlier combinations
  double t_SPA_Cutoff,                       // P-value cutoff for applying SPA correction
  double t_zeta,                             // SPA parameter for moment approximation
  double t_tol                               // Numerical tolerance for SPA convergence
) {
  if (ptr_gSPAGRMobj)
    delete ptr_gSPAGRMobj;

  ptr_gSPAGRMobj = new SPAGRM::SPAGRMClass(
    t_resid,                                 // Residual vector from null model
    t_resid_unrelated_outliers,              // Residuals for unrelated outlier subjects
    t_sum_R_nonOutlier,                      // Sum of residuals for non-outlier subjects
    t_R_GRM_R_nonOutlier,                    // Quadratic form R'*GRM*R for non-outliers
    t_R_GRM_R_TwoSubjOutlier,                // Quadratic form for two-subject outlier pairs
    t_R_GRM_R,                               // Full quadratic form R'*GRM*R
    t_MAF_interval,                          // Minor allele frequency intervals for binning
    t_TwoSubj_list,                          // List of two-subject outlier pairs information
    t_ThreeSubj_list,                        // List of three-subject outlier combinations
    t_SPA_Cutoff,                            // P-value cutoff for applying SPA correction
    t_zeta,                                  // SPA parameter for moment approximation
    t_tol                                    // Numerical tolerance for SPA convergence
  );
}

// Initialize SAGELD object for Scalable and Accurate Genome-wide Efficient Linear mixed model with Dosage.
// Handles gene-environment interactions with comprehensive outlier detection and SPA correction.
// [[Rcpp::export]]
void setSAGELDobjInCPP(
  std::string t_Method,                      // Analysis method ("SAIGE-Gene", "SAIGE-Gene+", etc.)
  arma::mat t_XTs,                           // Transpose of design matrix X
  arma::mat t_SS,                            // Matrix S'*S for efficient computation
  arma::mat t_AtS,                           // Matrix A'*S
  arma::mat t_Q,                             // Q matrix from QR decomposition
  arma::mat t_A21,                           // A21 matrix component
  arma::mat t_TTs,                           // Matrix T'*T
  arma::mat t_Tys,                           // Matrix T'*y
  arma::vec t_sol,                           // Solution vector from null model
  arma::vec t_blups,                         // Best linear unbiased predictors
  double t_sig,                              // Residual variance estimate
  arma::vec t_resid,                         // Main effect residuals
  arma::vec t_resid_G,                       // Genetic effect residuals
  arma::vec t_resid_GxE,                     // Gene-environment interaction residuals
  arma::vec t_resid_E,                       // Environmental effect residuals
  arma::vec t_resid_unrelated_outliers,      // Unrelated outlier residuals (main)
  arma::vec t_resid_unrelated_outliers_G,    // Unrelated outlier residuals (genetic)
  arma::vec t_resid_unrelated_outliers_GxE,  // Unrelated outlier residuals (GxE)
  double t_sum_R_nonOutlier,                 // Sum of residuals for non-outliers (main)
  double t_sum_R_nonOutlier_G,               // Sum of residuals for non-outliers (genetic)
  double t_sum_R_nonOutlier_GxE,             // Sum of residuals for non-outliers (GxE)
  double t_R_GRM_R,                          // Full quadratic form R'*GRM*R (main)
  double t_R_GRM_R_G,                        // Full quadratic form (genetic)
  double t_R_GRM_R_GxE,                      // Full quadratic form (GxE)
  double t_R_GRM_R_G_GxE,                    // Cross quadratic form (genetic-GxE)
  double t_R_GRM_R_E,                        // Full quadratic form (environmental)
  double t_R_GRM_R_nonOutlier,               // Non-outlier quadratic form (main)
  double t_R_GRM_R_nonOutlier_G,             // Non-outlier quadratic form (genetic)
  double t_R_GRM_R_nonOutlier_GxE,           // Non-outlier quadratic form (GxE)
  double t_R_GRM_R_nonOutlier_G_GxE,         // Non-outlier cross quadratic form
  double t_R_GRM_R_TwoSubjOutlier,           // Two-subject outlier quadratic form (main)
  double t_R_GRM_R_TwoSubjOutlier_G,         // Two-subject outlier quadratic form (genetic)
  double t_R_GRM_R_TwoSubjOutlier_GxE,       // Two-subject outlier quadratic form (GxE)
  double t_R_GRM_R_TwoSubjOutlier_G_GxE,     // Two-subject outlier cross quadratic form
  Rcpp::List t_TwoSubj_list,                 // Information for two-subject outlier pairs
  Rcpp::List t_ThreeSubj_list,               // Information for three-subject outlier combinations
  arma::vec t_MAF_interval,                  // Minor allele frequency intervals
  double t_zScoreE_cutoff,                   // Z-score cutoff for environmental effects
  double t_SPA_Cutoff,                       // P-value cutoff for SPA correction
  double t_zeta,                             // SPA parameter for moment approximation
  double t_tol                               // Numerical tolerance for convergence
) {
  if (ptr_gSAGELDobj)
    delete ptr_gSAGELDobj;

  ptr_gSAGELDobj = new SAGELD::SAGELDClass(
    t_Method,                                // Analysis method ("SAIGE-Gene", "SAIGE-Gene+", etc.)
    t_XTs,                                   // Transpose of design matrix X
    t_SS,                                    // Matrix S'*S for efficient computation
    t_AtS,                                   // Matrix A'*S
    t_Q,                                     // Q matrix from QR decomposition
    t_A21,                                   // A21 matrix component
    t_TTs,                                   // Matrix T'*T
    t_Tys,                                   // Matrix T'*y
    t_sol,                                   // Solution vector from null model
    t_blups,                                 // Best linear unbiased predictors
    t_sig,                                   // Residual variance estimate
    t_resid,                                 // Main effect residuals
    t_resid_G,                               // Genetic effect residuals
    t_resid_GxE,                             // Gene-environment interaction residuals
    t_resid_E,                               // Environmental effect residuals
    t_resid_unrelated_outliers,              // Unrelated outlier residuals (main)
    t_resid_unrelated_outliers_G,            // Unrelated outlier residuals (genetic)
    t_resid_unrelated_outliers_GxE,          // Unrelated outlier residuals (GxE)
    t_sum_R_nonOutlier,                      // Sum of residuals for non-outliers (main)
    t_sum_R_nonOutlier_G,                    // Sum of residuals for non-outliers (genetic)
    t_sum_R_nonOutlier_GxE,                  // Sum of residuals for non-outliers (GxE)
    t_R_GRM_R,                               // Full quadratic form R'*GRM*R (main)
    t_R_GRM_R_G,                             // Full quadratic form (genetic)
    t_R_GRM_R_GxE,                           // Full quadratic form (GxE)
    t_R_GRM_R_G_GxE,                         // Cross quadratic form (genetic-GxE)
    t_R_GRM_R_E,                             // Full quadratic form (environmental)
    t_R_GRM_R_nonOutlier,                    // Non-outlier quadratic form (main)
    t_R_GRM_R_nonOutlier_G,                  // Non-outlier quadratic form (genetic)
    t_R_GRM_R_nonOutlier_GxE,                // Non-outlier quadratic form (GxE)
    t_R_GRM_R_nonOutlier_G_GxE,              // Non-outlier cross quadratic form
    t_R_GRM_R_TwoSubjOutlier,                // Two-subject outlier quadratic form (main)
    t_R_GRM_R_TwoSubjOutlier_G,              // Two-subject outlier quadratic form (genetic)
    t_R_GRM_R_TwoSubjOutlier_GxE,            // Two-subject outlier quadratic form (GxE)
    t_R_GRM_R_TwoSubjOutlier_G_GxE,          // Two-subject outlier cross quadratic form
    t_TwoSubj_list,                          // Information for two-subject outlier pairs
    t_ThreeSubj_list,                        // Information for three-subject outlier combinations
    t_MAF_interval,                          // Minor allele frequency intervals
    t_zScoreE_cutoff,                        // Z-score cutoff for environmental effects
    t_SPA_Cutoff,                            // P-value cutoff for SPA correction
    t_zeta,                                  // SPA parameter for moment approximation
    t_tol                                    // Numerical tolerance for convergence
  );
}

// Initialize SPAmix object for saddle point approximation with population stratification correction.
// Uses principal components to account for population structure in mixed ancestry samples.
// [[Rcpp::export]]
void setSPAmixobjInCPP(
  arma::mat t_resid,                         // Matrix of residuals from null model
  arma::mat t_PCs,                           // Principal components matrix for population structure
  int t_N,                                   // Sample size
  double t_SPA_Cutoff,                       // P-value cutoff for applying SPA correction
  Rcpp::List t_outlierList                   // List containing outlier subject information
) {
  if (ptr_gSPAmixobj)
    delete ptr_gSPAmixobj;

  ptr_gSPAmixobj = new SPAmix::SPAmixClass(
    t_resid,                                 // Matrix of residuals from null model
    // t_XinvXX,                             // (Commented out) Inverse design matrix
    // t_tX,                                 // (Commented out) Transpose design matrix
    t_PCs,                                   // Principal components matrix for population structure
    t_N,                                     // Sample size
    t_SPA_Cutoff,                            // P-value cutoff for applying SPA correction
    t_outlierList                            // List containing outlier subject information
  );
}

// Initialize SPACox object for saddle point approximation in Cox proportional hazards models.
// Handles time-to-event analysis with efficient variance estimation and SPA correction.
// [[Rcpp::export]]
void setSPACoxobjInCPP(
  arma::mat t_cumul,                         // Cumulative hazard matrix
  arma::vec t_mresid,                        // Martingale residuals from Cox model
  arma::mat t_XinvXX,                        // (X'X)^(-1) matrix for variance calculation
  arma::mat t_tX,                            // Transpose of design matrix X
  int t_N,                                   // Sample size
  double t_pVal_covaAdj_Cutoff,              // P-value cutoff for covariate adjustment
  double t_SPA_Cutoff                        // P-value cutoff for applying SPA correction
) {
  if (ptr_gSPACoxobj)
    delete ptr_gSPACoxobj;

  ptr_gSPACoxobj = new SPACox::SPACoxClass(
    t_cumul,                                 // Cumulative hazard matrix
    t_mresid,                                // Martingale residuals from Cox model
    t_XinvXX,                                // (X'X)^(-1) matrix for variance calculation
    t_tX,                                    // Transpose of design matrix X
    t_N,                                     // Sample size
    t_pVal_covaAdj_Cutoff,                   // P-value cutoff for covariate adjustment
    t_SPA_Cutoff                             // P-value cutoff for applying SPA correction
  );
}

// Initialize WtCoxG object for weighted Cox regression analysis with genetic variants.
// Handles survival analysis with saddle point approximation and external reference corrections.
// [[Rcpp::export]]
void setWtCoxGobjInCPP(
  Rcpp::DataFrame t_mergeGenoInfo,           // Merged genotype information for batch effect testing
  arma::vec t_mresid,                        // Martingale residuals from Cox model
  arma::vec t_weight,                        // Weight vector for analysis
  double t_cutoff,                           // batch effect p-value cutoff for association testing
  double t_SPA_Cutoff                        // P-value cutoff for applying SPA correction
) {
  if (ptr_gWtCoxGobj)
    delete ptr_gWtCoxGobj;

  ptr_gWtCoxGobj = new WtCoxG::WtCoxGClass(
    t_mergeGenoInfo,                         // Merged genotype information for batch effect testing
    t_mresid,                                // Martingale residuals from Cox model (R vector)
    t_weight,                                // Weight vector for analysis (w vector)
    t_cutoff,                                // batch effect p-value cutoff for association testing
    t_SPA_Cutoff                             // P-value cutoff for applying SPA correction
  );
}
