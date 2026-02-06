#' Instruction of SPAmixLocalPlus method
#'
#' SPAmixLocalPlus performs local ancestry-specific association testing by accounting 
#' for local ancestry patterns using ancestry-specific haplotype counts and pre-computed 
#' phi (intuitively kinship coefficients stratified by local ancestry) matrices.
#' 
#' @return NULL
#'
#' @examples
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' GRM <- data.table::fread(SparseGRMFile)
#'
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile)
#' subjData <- PhenoData$IID
#'
#' extdata_dir <- system.file("extdata", package = "GRAB")
#' dosagePrefix <- paste0(extdata_dir, "/simuAncestry")
#' 
#' temp_dir <- tempdir()
#' phiOutputPrefix <- paste0(temp_dir, "/phi_")
#' outPrefix <- paste0(temp_dir, "/out_")
#'
#' # Step0: Estimate phi for ancestries 1 and 2
#' SPAmixLocalPlus.EstimatePhi(
#' GRM = GRM,
#'   dosagePrefix = dosagePrefix,
#'   haploPrefix = NULL,  # Uses dosagePrefix if NULL
#'   ancIdx = c(1, 2),
#'   SampleIDs = subjData,
#'   phiOutputPrefix = phiOutputPrefix,
#'   MAF_cutoff = 0.01
#' )
#'
#' # Step1: Fit null model to get residuals
#' residuals <- survival::coxph(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData
#' )$residuals
#'
#' # Step2: Run SPAmixPlus tests for ancestries 1 and 2
#' result_lst <- SPAmixLocalPlus.Marker(
#'   resid = residuals,
#'   subjData = subjData,
#'   dosagePrefix = dosagePrefix,
#'   haploPrefix = NULL,  # Uses dosagePrefix if NULL
#'   phiPrefix = phiOutputPrefix,
#'   ancIdx = c(1, 2),
#'   outPrefix = outPrefix
#' )
#' 
#' # View results
#' OutputFile1 <- paste0(outPrefix, "Ancestry1.txt")
#' head(data.table::fread(OutputFile1))
#'
#' @seealso
#' \code{\link{SPAmixLocalPlus.EstimatePhi}} for phi matrix estimation,
#' \code{\link{SPAmixLocalPlus.Marker}} for the main analysis function
#' 
#' @references
#' Yuzhuo Ma (2026). \url{https://github.com/YuzhuoMa97/SPAmixPlus}
#'
GRAB.SPAmixLocalPlus <- function() {
  .message("?GRAB.SPAmixLocalPlus for instructions")
}

#' Introduction of SPAmixLocalPlus.EstimatePhi
#'
#' Estimates ancestry-specific kinship coefficients (phi matrices) for local ancestry
#' association testing. Computes pairwise genetic correlations stratified by haplotype
#' count configurations (scenarios A, B, C, D) using dosage and haplotype count data.
#' 
#' @return NULL
#'
#' @param GRM Three columns: ID1, ID2, Value.
#'   Contains pairwise kinship coefficients for related individuals. Typically from
#'   \code{getSparseGRM()} or similar GRM estimation.
#' @param dosagePrefix Character (required). Prefix path for dosage files. Function expects
#'   files named \code{paste0(dosagePrefix, ancIdx, "_Dosage.txt.gz")} for each ancestry.
#' @param haploPrefix Character (optional, default = \code{NULL}). Prefix path for haplotype
#'   count files. If \code{NULL}, uses \code{dosagePrefix}. Function expects files named
#'   \code{paste0(haploPrefix, ancIdx, "_HapCount.txt.gz")} for each ancestry.
#' @param ancIdx Integer vector (required). Ancestry indices to process (e.g., \code{c(1, 2)}).
#' @param SampleIDs Character vector (required). Sample IDs matching those in dosage/hapcount files.
#'   Must match order in \code{GRM}.
#' @param phiOutputPrefix Character (required). Prefix path for phi output files. Function will create
#'   files named \code{paste0(phiOutputPrefix, "Ancestry", ancIdx, "_scenario", c("A","B","C","D"), ".txt")}.
#'   Example: if \code{phiOutputPrefix = "./phi/result_"}, creates \code{./phi/result_Ancestry1_scenarioA.txt}, etc.
#' @param Scenarios Character vector (optional, default = \code{c("A", "B", "C", "D")}).
#'   Phi scenarios to compute. Each scenario corresponds to haplotype count configurations:
#'   \itemize{
#'     \item A: Both individuals have 2 haplotypes from this ancestry
#'     \item B: Individual i has 2, individual j has 1 haplotype
#'     \item C: Individual i has 1, individual j has 2 haplotypes
#'     \item D: Both individuals have 1 haplotype from this ancestry
#'   }
#' @param Threshold Numeric (optional, default = 0.0). Phi threshold for filtering output.
#'   Only pairs with phi >= Threshold are saved.
#' @param MAF_cutoff Numeric (optional, default = 0.01). Minor Allele Frequency cutoff.
#'   SNPs with MAF < MAF_cutoff or MAF > (1 - MAF_cutoff) are excluded from phi estimation.
#' 
#' \deqn{phi_{ij} = \frac{1}{M} \sum_{s=1}^{M} \frac{(g_{is} - h_{is}q_s)(g_{js} - h_{js}q_s)}{h_{is} h_{js} q_s(1-q_s)}}
#' 
#' where:
#' \itemize{
#'   \item \eqn{g_{is}} = dosage for individual i at SNP s
#'   \item \eqn{h_{is}} = haplotype count for individual i at SNP s
#'   \item \eqn{q_s} = global allele frequency at SNP s
#'   \item M = number of SNPs passing filters
#' }
#'
#' @return NULL Results are saved with filenames:
#'   \code{{phiOutputPrefix}Ancestry{ancIdx}_scenario{A/B/C/D}.txt}
#'   
#'   Each file is tab-delimited with columns:
#'   \itemize{
#'     \item i: Sample ID for individual i
#'     \item j: Sample ID for individual j
#'     \item phi_value: Estimated phi coefficient
#'   }
#'   
#'   Files contain bidirectional pairs (i,j) and (j,i) for all pairs in GRM.
#'
SPAmixLocalPlus.EstimatePhi <- function(
  GRM,
  dosagePrefix,
  haploPrefix = NULL,
  ancIdx,
  SampleIDs,
  phiOutputPrefix,
  Scenarios = c("A", "B", "C", "D"),
  Threshold = 0.0,
  MAF_cutoff = 0.01
) {
  
  # Setup output directory from prefix
  OutputDir <- dirname(phiOutputPrefix)
  if (OutputDir != "." && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
  
  # Process GRM
  if (is.character(GRM)) {
    if (!file.exists(GRM)) stop("GRM file not found: ", GRM)
    grm_dt <- data.table::fread(GRM)
  } else {
    grm_dt <- data.table::as.data.table(GRM)
  }
  
  # Standardize GRM column names
  if (ncol(grm_dt) < 3) stop("GRM must have at least 3 columns (ID1, ID2, Value).")
  if (!all(c("ID1", "ID2") %in% names(grm_dt))) {
    colnames(grm_dt)[1:3] <- c("ID1", "ID2", "Value")
  }
  
  # Map sample IDs to 0-based indices for C++
  grm_dt$idx1 <- match(grm_dt$ID1, SampleIDs) - 1
  grm_dt$idx2 <- match(grm_dt$ID2, SampleIDs) - 1
  
  valid_grm <- grm_dt[!is.na(grm_dt$idx1) & !is.na(grm_dt$idx2), ]
  if (nrow(valid_grm) < nrow(grm_dt)) {
    warning("Some GRM pairs contain IDs not in SampleIDs list.")
  }
  
  if (is.null(haploPrefix)) {
    haploPrefix <- dosagePrefix
  }
  
  # Process each ancestry
  for (anc_id in ancIdx) {
    dosage_file <- paste0(dosagePrefix, anc_id, "_Dosage.txt.gz")
    haplo_file <- paste0(haploPrefix, anc_id, "_HapCount.txt.gz")
    
    if (!file.exists(dosage_file)) {
      warning("Dosage file not found, skipping ancestry ", anc_id, ": ", dosage_file)
      next
    }
    if (!file.exists(haplo_file)) {
      warning("Haplotype file not found, skipping ancestry ", anc_id, ": ", haplo_file)
      next
    }
    
    .message("Processing Ancestry %s...", anc_id)
    
    # Read dosage and haplotype count files
    .message("Reading dosage file: %s", dosage_file)
    dos_dt <- data.table::fread(dosage_file, header = TRUE)
    dos_mat <- as.matrix(dos_dt[, -c(1:5)])  # Remove first 5 columns (CHROM, POS, ID, REF, ALT)
    rm(dos_dt)
    .message("Reading dosage file completed: %d SNPs x %d samples", nrow(dos_mat), ncol(dos_mat))
    
    .message("Reading haplotype count file: %s", haplo_file)
    hap_dt <- data.table::fread(haplo_file, header = TRUE)
    hap_mat <- as.matrix(hap_dt[, -c(1:5)])  # Remove first 5 columns
    rm(hap_dt)
    .message("Reading haplotype count file completed")
    
    # Process each scenario
    for (scenario in Scenarios) {
      .message("Analyzing Scenario %s...", scenario)
      
      # Call C++ function
      res <- SPAmixLocalPlus_computePhiInCPP(
        hapcount_matrix = hap_mat,
        dosage_matrix = dos_mat,
        pair_idx1 = valid_grm$idx1,
        pair_idx2 = valid_grm$idx2,
        scenario = scenario,
        phi_threshold = Threshold,
        maf_cutoff = MAF_cutoff
      )
      
      # Compute final phi values
      Phi_Values <- res$ratio_sums / pmax(res$valid_counts, 1)
      Phi_Values[res$valid_counts == 0] <- NA
      
      # C++ function already returns bidirectional pairs, use them directly
      final_pair_idx1 <- res$directed_pair_idx1
      final_pair_idx2 <- res$directed_pair_idx2
      
      # Convert indices to sample IDs (C++ uses 0-based indexing)
      final_ids_i <- SampleIDs[final_pair_idx1 + 1]
      final_ids_j <- SampleIDs[final_pair_idx2 + 1]
      
      res_df <- data.table::data.table(
        i = final_ids_i,
        j = final_ids_j,
        phi_value = Phi_Values
      )
      
      # Filter out NA and apply threshold
      res_df <- res_df[!is.na(res_df$phi_value), ]
      if (Threshold > 0) res_df <- res_df[res_df$phi_value >= Threshold, ]
      
      # Save output with compatible naming: {phiOutputPrefix}Ancestry{ancIdx}_scenario{A/B/C/D}.txt
      out_path <- paste0(phiOutputPrefix, "Ancestry", anc_id, "_scenario", scenario, ".txt")
      
      data.table::fwrite(res_df, out_path, sep = "\t", quote = FALSE)
      .message("Scenario %s saved: %d pairs", scenario, nrow(res_df))
    }
  }
  
  .message("Phi estimation completed.")
  invisible(NULL)
}

#' Introduction of SPAmixLocalPlus.Marker
#'
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
#' @param phiPrefix Character (optional, default = \code{NULL}). Prefix path for phi matrix files.
#'   If \code{NULL}, uses \code{dosagePrefix}.
#'   Function expects files named \code{paste0(phiPrefix, "Ancestry", ancIdx, "_scenario", c("A","B","C","D"), ".txt")}
#'   for each ancestry. Should match the \code{phiOutputPrefix} used in \code{SPAmixLocalPlus.EstimatePhi()}.
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
#' @return NULL Results are saved with filenames:
#'   \code{paste0(outPrefix, "Ancestry", ancIdx, ".txt")} for each ancestry.
#'
#' @seealso
#' \code{\link{GRAB.SPAmixLocalPlus}} for usage examples,
#' \code{\link{SPAmixLocalPlus.OneAnc}} for single-ancestry processing
#'
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

  for (i in ancIdx) {
    SPAmixLocalPlus.OneAnc(
      dosage_file = paste0(dosagePrefix, i, "_Dosage.txt.gz"),
      haplo_file  = paste0(haploPrefix, i, "_HapCount.txt.gz"),
      phi_anc_pre = paste0(phiPrefix, "Ancestry", i, "_scenario"),
      output_file = paste0(outPrefix, "Ancestry", i, ".txt")
    )
  }

  return(NULL)
}


#' Introdution of SPAmixLocalPlus.OneAnc
#' 
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
#' 
#' @return NULL Results are also written to \code{output_file}
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
  matched_sample_ids <- file_sample_ids[match_info$file_indices + 1]
  
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
      dt = dt[!is.na(dt$i) & !is.na(dt$j), ] # Keep only pairs in our sample
      
      if(nrow(dt) == 0) return(matrix(0, 0, 3))
      
      # Return matrix (i, j, value)
      as.matrix(dt[, c("i", "j", "val")])
  }
  
  phi_A_mat <- convert_phi(paste0(phi_anc_pre, "A.txt"), matched_sample_ids)
  phi_B_mat <- convert_phi(paste0(phi_anc_pre, "B.txt"), matched_sample_ids)
  phi_C_mat <- convert_phi(paste0(phi_anc_pre, "C.txt"), matched_sample_ids)
  phi_D_mat <- convert_phi(paste0(phi_anc_pre, "D.txt"), matched_sample_ids)
  
  # 4. Run Analysis (using global state set in SPAmixPlusLocal_setupInCPP)
  SPAmixPlusLocal_streamInCPP(
    geno_file = dosage_file,
    haplo_file = haplo_file, 
    output_file = output_file,
    file_match_idx = match_info$file_indices,
    phi_A_mat = phi_A_mat,
    phi_B_mat = phi_B_mat,
    phi_C_mat = phi_C_mat,
    phi_D_mat = phi_D_mat
  )
  
  return(NULL)
}
