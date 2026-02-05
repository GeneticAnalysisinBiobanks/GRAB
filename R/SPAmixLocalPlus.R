#' Local Ancestry-Specific Association Testing Using SPAmixLocalPlus
#'
#' @description
#' Main entry point for SPAmixLocalPlus method. This function provides instructions
#' for performing local ancestry-specific association testing in admixed populations.
#' Unlike global ancestry methods, SPAmixLocalPlus accounts for local ancestry patterns
#' using ancestry-specific haplotype counts and pre-computed phi matrices (pairwise
#' genetic correlation coefficients).
#'
#' @details
#' SPAmixLocalPlus extends SPAmixPlus to handle local ancestry-specific genetic associations.
#' The workflow involves:
#' \enumerate{
#'   \item Prepare ancestry-specific dosage files (*.txt.gz) with ALT allele dosages
#'   \item Prepare ancestry-specific haplotype count files (*.txt.gz)
#'   \item Pre-compute phi matrices for scenarios A, B, C, D per ancestry
#'   \item Fit null model using residuals from survival or regression model
#'   \item Run \code{SPAmixLocalPlus.Marker()} to analyze all ancestries
#' }
#'
#' @return No return value. Prints usage instructions to console.
#'
#' @examples
#' \dontrun{
#' # Load phenotype data
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile)
#' subjData <- PhenoData$IID
#'
#' # Fit null model and extract residuals
#' residuals <- survival::coxph(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData
#' )$residuals
#'
#' # Set file paths
#' dosagePrefix <- system.file("extdata", "simuAncestry", package = "GRAB")
#' # Expects: simuAncestry1_Dosage.txt.gz, simuAncestry1_HapCount.txt.gz, etc.
#' # and: simuAncestryAncestry1_scenarioA.txt, etc.
#'
#' outPrefix <- file.path(tempdir(), "resultSPAmixLocalPlus_Ancestry")
#'
#' # Run analysis for ancestries 1 and 2
#' result_lst <- SPAmixLocalPlus.Marker(
#'   resid = residuals,
#'   subjData = subjData,
#'   dosagePrefix = dosagePrefix,
#'   haploPrefix = NULL,  # Uses dosagePrefix if NULL
#'   phiPrefix = NULL,    # Uses dosagePrefix if NULL
#'   ancIdx = c(1, 2),
#'   outPrefix = outPrefix
#' )
#'
#' # View results for Ancestry 1
#' OutputFile1 <- paste0(outPrefix, "Ancestry1.txt")
#' head(data.table::fread(OutputFile1))
#' }
#'
#' @seealso
#' \code{\link{SPAmixLocalPlus.Marker}} for the main analysis function
#'
#' @export
GRAB.SPAmixLocalPlus <- function() {
  .message("?GRAB.SPAmixLocalPlus for instructions")
}

#' Run SPAmixLocalPlus Analysis for Multiple Ancestries
#'
#' @description
#' Main worker function to perform local ancestry-specific association testing across
#' multiple ancestries. Sets global state once, then processes each ancestry sequentially.
#'
#' @param resid Numeric vector (required). Residuals from null model fit (e.g., from
#'   \code{survival::coxph()$residuals} or \code{lm()$residuals}). Must match order of \code{subjData}.
#' @param subjData Character vector (required). Sample IDs corresponding to \code{resid}.
#'   Used to match samples between phenotype data and genotype files.
#' @param dosagePrefix Character (required). Prefix path for dosage files. Function expects
#'   files named \code{paste0(dosagePrefix, ancIdx, "_Dosage.txt.gz")} for each ancestry.
#'   Example: if \code{dosagePrefix = "data/Anc"} and \code{ancIdx = c(1,2)}, expects
#'   \code{data/Anc1_Dosage.txt.gz} and \code{data/Anc2_Dosage.txt.gz}.
#' @param haploPrefix Character (optional, default = \code{NULL}). Prefix path for haplotype
#'   count files. If \code{NULL}, uses \code{dosagePrefix}. Function expects files named
#'   \code{paste0(haploPrefix, ancIdx, "_HapCount.txt.gz")} for each ancestry.
#' @param phiPrefix Character (optional, default = \code{NULL}). Prefix path for phi matrix
#'   files. If \code{NULL}, uses \code{dosagePrefix}. Function expects files named
#'   \code{paste0(phiPrefix, "Ancestry", ancIdx, "_scenario", c("A","B","C","D"), ".txt")}
#'   for each ancestry.
#' @param ancIdx Integer vector (required). Ancestry indices to analyze (e.g., \code{c(1, 2, 3)}).
#'   Each index corresponds to one local ancestry population.
#' @param outPrefix Character (required). Prefix path for output files. Results saved as
#'   \code{paste0(outPrefix, "Ancestry", ancIdx, ".txt")} for each ancestry.
#' @param outLierList Character vector (optional, default = \code{NULL}). Sample IDs to
#'   exclude from analysis as outliers.
#' @param save_interval Integer (optional, default = 100). Number of SNPs processed before
#'   intermediate results are saved to disk. Helps prevent data loss in long-running analyses.
#' @param MAF_cutoff Numeric (optional, default = 0.00001). Minor Allele Frequency cutoff.
#'   SNPs with MAF below this threshold are excluded from analysis.
#' @param MAC_cutoff Numeric (optional, default = 1). Minor Allele Count cutoff. SNPs with
#'   MAC below this threshold are excluded from analysis.
#' @param cutoff Numeric (optional, default = 2.0). Saddlepoint approximation cutoff for
#'   p-value calculation. Uses normal approximation when |test statistic| < cutoff.
#' @param verbose Logical (optional, default = \code{TRUE}). If \code{TRUE}, prints detailed
#'   progress information during analysis.
#'
#' @return List with one element per ancestry (indexed by \code{ancIdx}). Each element
#'   contains analysis results from \code{SPAmixLocalPlus.OneAnc()}. Results are also
#'   written to output files.
#'
#' @seealso
#' \code{\link{GRAB.SPAmixLocalPlus}} for usage examples,
#' \code{\link{SPAmixLocalPlus.OneAnc}} for single-ancestry processing
#'
#' @export
SPAmixLocalPlus.Marker <- function(
  resid,
  subjData,
  dosagePrefix,
  haploPrefix = NULL,
  phiPrefix = NULL,
  ancIdx,
  outPrefix,
  outLierList = NULL,
  save_interval = 100,
  MAF_cutoff = 0.00001,
  MAC_cutoff = 1,
  cutoff = 2.0,
  verbose = TRUE
) {

  resid <- as.numeric(resid)

  # Set global state ONCE for all ancestries
  SPAmixPlusLocal_setupInCPP(
    resid = resid,
    subjData = subjData,
    outLierList = outLierList,
    save_interval = save_interval,
    MAF_cutoff = MAF_cutoff,
    MAC_cutoff = MAC_cutoff,
    cutoff = cutoff,
    verbose = verbose
  )

  if (is.null(haploPrefix)) {
    haploPrefix <- dosagePrefix
  }

  if (is.null(phiPrefix)) {
    phiPrefix <- dosagePrefix
  }

  result_lst <- vector("list", length(ancIdx))
  for (i in ancIdx) {
    result_lst[[i]] <- SPAmixLocalPlus.OneAnc(
      dosage_file = paste0(dosagePrefix, i, "_Dosage.txt.gz"),
      haplo_file  = paste0(haploPrefix, i, "_HapCount.txt.gz"),
      phi_anc_pre = paste0(phiPrefix, "Ancestry", i, "_scenario"),
      output_file = paste0(outPrefix, "Ancestry", i, ".txt")
    )
  }

  return(result_lst)
}


# ============================== SPAmixLocalPlus: Step2 =================================

#' Process One Ancestry for SPAmixLocalPlus Analysis
#' 
#' @description
#' Worker function to perform local ancestry-specific association testing for a single
#' ancestry. Uses global state set by \code{SPAmixLocalPlus.Marker()} for residuals,
#' sample IDs, and analysis parameters. Loads ancestry-specific dosage, haplotype counts,
#' and phi matrices to compute association statistics.
#' 
#' @param dosage_file Character (required). Full path to ancestry-specific genotype dosage
#'   file (.txt.gz). Contains dosage values (0-2) for ALT allele from specific ancestry.
#'   Format: tab-delimited, SNPs in rows, samples in columns, with SNP ID in first column
#'   and sample IDs in header row.
#' @param haplo_file Character (required). Full path to ancestry-specific haplotype count
#'   file (.txt.gz). Contains number of haplotypes (0, 1, or 2) from specific ancestry at
#'   each locus. Same format as dosage file.
#' @param phi_anc_pre Character (required). Prefix path to phi files for this ancestry
#'   (without scenario suffix). Function expects files: \code{paste0(phi_anc_pre, c("A","B","C","D"), ".txt")}.
#'   Example: if \code{phi_anc_pre = "phi/ancestry1_scenario"}, expects
#'   \code{phi/ancestry1_scenarioA.txt}, \code{phi/ancestry1_scenarioB.txt}, etc.
#'   Phi files should be tab-delimited with columns: IID1, IID2, Phi_value.
#' @param output_file Character (required). Full path to output result file where
#'   association test results will be saved (.txt format).
#' 
#' @details 
#' The \code{SPAmixLocalPlus.Marker} function performs association testing that accounts for 
#' local ancestry patterns. Unlike global ancestry methods, it uses:
#' 
#' \itemize{
#'   \item \strong{Dosage files}: Ancestry-specific allele dosages
#'   \item \strong{Haplotype counts}: Number of haplotypes from each ancestry at each locus
#'   \item \strong{Phi matrices}: Pairwise genetic correlations for different haplotype scenarios
#' }
#' 
#' ## Phi Files and Scenarios
#' 
#' Phi files contain pairwise genetic correlation coefficients specific to local ancestry,
#' stratified by haplotype count configurations:
#' 
#' \itemize{
#'   \item \strong{Scenario A}: Both individuals have 2 haplotypes from this ancestry
#'   \item \strong{Scenario B}: Individual i has 2, individual j has 1 haplotype
#'   \item \strong{Scenario C}: Individual i has 1, individual j has 2 haplotypes
#'   \item \strong{Scenario D}: Both individuals have 1 haplotype from this ancestry
#' }
#' 
#' ## File Format
#' 
#' Input files should be tab-separated with samples in columns and SNPs in rows:
#' - First column: SNP ID
#' - Header row: Sample IDs
#' - Data: Numeric values (dosages or haplotype counts)
#' 
#' @return List containing analysis results. Results are also written to \code{output_file}
#'   as tab-delimited text with columns:
#'   \itemize{
#'     \item Marker ID
#'     \item Chromosome
#'     \item Position
#'     \item Alleles
#'     \item MAF (Minor Allele Frequency)
#'     \item MAC (Minor Allele Count)
#'     \item P-value
#'     \item Other test statistics
#'   }
#' 
#' @note
#' This function is typically called by \code{SPAmixLocalPlus.Marker()} rather than directly.
#' It relies on global state variables set by \code{SPAmixPlusLocal_setupInCPP()}, including
#' residuals, sample IDs, outlier list, and analysis parameters.
#' 
#' @seealso
#' \code{\link{SPAmixLocalPlus.Marker}} for multi-ancestry analysis,
#' \code{\link{GRAB.SPAmixLocalPlus}} for usage examples
#' 
#' @export

SPAmixLocalPlus.OneAnc = function(
  dosage_file,
  haplo_file, 
  phi_anc_pre,
  output_file
) {
  
  # 1. Validation - Check file existence
  if (!file.exists(dosage_file)) stop("Dosage file not found: ", dosage_file)
  if (!file.exists(haplo_file)) stop("Haplotype file not found: ", haplo_file)
  
  # Check for phi files (scenarios A, B, C, D)
  phi_scenarios <- c("A", "B", "C", "D")
  for (scenario in phi_scenarios) {
    phi_path <- paste0(phi_anc_pre, scenario, ".txt")
    if (!file.exists(phi_path)) {
      stop("Phi file not found: ", phi_path, 
           "\nExpected files: ", phi_anc_pre, "A/B/C/D.txt")
    }
  }
  
  # 2. Match samples and get indices from C++ global state
  file_sample_ids <- read_ukb_sample_ids_cpp(dosage_file)
  match_info <- getSampleMatchIndices_cpp(file_sample_ids)
  
  # 3. Load Phi Matrices
  # Helper to load and format phi matrix using matched sample indices
  # Get matched sample IDs for phi mapping
  matched_sample_ids <- match_info$matched_sample_ids
  
  convert_phi <- function(phi_path, sample_ids) {
      # Read phi file
      p_data <- data.table::fread(phi_path, header = TRUE)
      
      if(nrow(p_data) == 0) return(matrix(0, 0, 3))
      
      # Assume columns: IID1, IID2, Phi_value (or similar names)
      colnames(p_data) <- c("i_id", "j_id", "phi_value")
      
      i_ids = as.character(p_data$i_id)
      j_ids = as.character(p_data$j_id)
      vals = as.numeric(p_data$phi_value)
      
      # Create mapping from sample IDs to indices
      sample_map_idx <- setNames(seq_along(sample_ids), sample_ids)
      
      # Map to internal indices (1-based for R)
      idx_i = sample_map_idx[i_ids]
      idx_j = sample_map_idx[j_ids]
      
      dt = data.table::data.table(i = idx_i, j = idx_j, val = vals)
      dt = dt[!is.na(i) & !is.na(j)] # Keep only pairs in our sample
      
      if(nrow(dt) == 0) return(matrix(0, 0, 3))
      
      # Return matrix (i, j, value)
      as.matrix(dt[, .(i, j, val)])
  }
  
  phi_A_mat <- convert_phi(paste0(phi_anc_pre, "A.txt"), matched_sample_ids)
  phi_B_mat <- convert_phi(paste0(phi_anc_pre, "B.txt"), matched_sample_ids)
  phi_C_mat <- convert_phi(paste0(phi_anc_pre, "C.txt"), matched_sample_ids)
  phi_D_mat <- convert_phi(paste0(phi_anc_pre, "D.txt"), matched_sample_ids)
  
  # 4. Run Analysis (using global state set in SPAmixPlusLocal_setupInCPP)
  result <- SPAmixPlusLocal_streamInCPP(
    geno_file = dosage_file,
    haplo_file = haplo_file, 
    output_file = output_file,
    file_match_idx = match_info$file_indices,
    phi_A_mat = phi_A_mat,
    phi_B_mat = phi_B_mat,
    phi_C_mat = phi_C_mat,
    phi_D_mat = phi_D_mat
  )
  
  return(result)
}
