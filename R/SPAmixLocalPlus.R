

# ============================== SPAmixLocalPlus: Step2 =================================

#' SPAmixLocalPlus Marker Analysis for Local Ancestry-Specific Association Testing
#' 
#' @description
#' Performs local ancestry-specific association testing using pre-computed phi matrices 
#' (pairwise genetic correlation coefficients). This method accounts for local ancestry 
#' patterns in admixed populations by using ancestry-specific haplotype counts and 
#' phi coefficients for different haplotype configuration scenarios.
#' 
#' @param objNull Object of class `SPAmixPlus_NULL_Model` returned by `fitNullModel.SPAmixPlus()`.
#'   Must contain residuals, sample IDs, outlier information, and other null model parameters.
#' @param dosage_file Character. Path to ancestry-specific genotype dosage file (.txt.gz).
#'   This file contains dosage values (0-2) for the ALT allele from a specific ancestry population.
#' @param haplo_file Character. Path to ancestry-specific haplotype count file (.txt.gz).
#'   This file contains the number of haplotypes (0, 1, or 2) from the specific ancestry at each locus.
#' @param phi_file Character. Prefix path to phi files for this ancestry (without scenario suffix).
#'   The function will look for files: phi_file_scenarioA.txt, phi_file_scenarioB.txt, etc.
#'   If phi_file = "phi/ancestry1", it will load "phi/ancestry1_scenarioA.txt", etc.
#' @param output_file Character. Path to output result file where association test results will be saved.
#' @param pheno_idx Integer. Phenotype index to analyze (default 1). Used when objNull contains multiple phenotypes.
#' @param save_interval Integer. Number of SNPs processed before intermediate results are saved (default 100).
#'   Helps prevent data loss for long-running analyses.
#' @param MAF_cutoff Numeric. Minor Allele Frequency cutoff (default 0.00001). SNPs with MAF below this are excluded.
#' @param MAC_cutoff Numeric. Minor Allele Count cutoff (default 1). SNPs with MAC below this are excluded.
#' @param verbose Logical. If TRUE, output detailed progress information (default TRUE).
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
#' @return No explicit return value. Results are written to \code{output_file} containing:
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
#' @examples
#' \dontrun{
#' # First, fit null model (same as SPAmixPlus)
#' objNull <- fitNullModel.SPAmixPlus(
#'   formula = pheno ~ age + sex,
#'   data = phenoData,
#'   SparseGRMFile = "grm/sparse_grm.txt",
#'   subjData = subj_ids
#' )
#' 
#' # Run local ancestry analysis for Ancestry 1
#' SPAmixLocalPlus.Marker(
#'   objNull = objNull,
#'   dosage_file = "data/Ancestry1_Dosage.txt.gz",
#'   haplo_file = "data/Ancestry1_HapCount.txt.gz",
#'   phi_file = "phi/ancestry1",  # Will look for ancestry1_scenarioA.txt, etc.
#'   output_file = "results/Ancestry1_results.txt"
#' )
#' 
#' # Run for Ancestry 2
#' SPAmixLocalPlus.Marker(
#'   objNull = objNull,
#'   dosage_file = "data/Ancestry2_Dosage.txt.gz",
#'   haplo_file = "data/Ancestry2_HapCount.txt.gz",
#'   phi_file = "phi/ancestry2",
#'   output_file = "results/Ancestry2_results.txt"
#' )
#' }
#' 
#' @seealso 
#' \code{\link{fitNullModel.SPAmixPlus}} for fitting the null model
#' 
#' @export
SPAmixLocalPlus.Marker = function(objNull,
                                   dosage_file,
                                   haplo_file, 
                                   phi_file,
                                   output_file,
                                   pheno_idx = 1,
                                   save_interval = 100,
                                   MAF_cutoff = 0.00001,
                                   MAC_cutoff = 1,
                                   verbose = TRUE) {
  
  # 1. Validation
  if (verbose) cat("=== SPAmixLocalPlus.Marker ===\n")
  if (!file.exists(dosage_file)) stop("Dosage file not found: ", dosage_file)
  if (!file.exists(haplo_file)) stop("Haplotype file not found: ", haplo_file)
  
  # Check for phi files (scenarios A, B, C, D)
  phi_scenarios <- c("A", "B", "C", "D")
  phi_file_pattern <- paste0(phi_file, "_scenario")
  for (scenario in phi_scenarios) {
    phi_path <- paste0(phi_file_pattern, scenario, ".txt")
    if (!file.exists(phi_path)) {
      stop("Phi file not found: ", phi_path, 
           "\nExpected files: ", phi_file, "_scenarioA/B/C/D.txt")
    }
  }
  
  # Check class of objNull
  # Note: The user might call it standard list or class. 
  # We check for basic required components if class check fails or behaves weirdly.
  if (!inherits(objNull, "SPAmixPlus_NULL_Model") && !inherits(objNull, "list")) {
      warning("objNull is not of class 'SPAmixPlus_NULL_Model' or 'list'. Proceeding with caution.")
  }

  if (verbose) {
    cat("Dosage File:", dosage_file, "\n")
    cat("Haplotype File:", haplo_file, "\n")
    cat("Phi File Prefix:", phi_file, "\n")
    cat("Output File:", output_file, "\n")
    cat("MAF threshold:", MAF_cutoff, "| MAC threshold:", MAC_cutoff, "\n")
  }

  # 2. Extract Null Model Info
  # Handle cases where nPheno might not be present (legacy objects?)
  nPheno = if(!is.null(objNull$nPheno)) objNull$nPheno else 1
  if (pheno_idx > nPheno || pheno_idx < 1) stop("pheno_idx out of range.")
  
  # Extract Residuals
  # objNull$resid is usually a matrix (N x nPheno)
  if(is.null(objNull$resid)) stop("objNull$resid is missing.")
  R = if(is.matrix(objNull$resid)) as.numeric(objNull$resid[, pheno_idx]) else as.numeric(objNull$resid)
  
  # Extract Sample IDs
  if (is.list(objNull$subjData) && "SubjID" %in% names(objNull$subjData)) {
      objNull_sample_ids <- as.character(objNull$subjData$SubjID)
  } else if(!is.null(objNull$subjData)) {
      objNull_sample_ids <- as.character(objNull$subjData)
  } else {
      # Fallback: try to find sample IDs from resid names?
      if(!is.null(rownames(objNull$resid))) {
          objNull_sample_ids = rownames(objNull$resid)
      } else {
          stop("Could not find Sample IDs in objNull$subjData or rownames(objNull$resid).")
      }
  }
  
  # Extract Outliers (0-based indices from C++)
  # objNull$outLierList is a list of lists if nPheno > 1, or just a list?
  # Standard SPAmixPlus usually lists of lists.
  outlier_list <- if(is.null(objNull$outLierList[[pheno_idx]])) objNull$outLierList else objNull$outLierList[[pheno_idx]]
  
  if(is.null(outlier_list$posOutlier)) {
      # Maybe structure is different
      posOutlier = integer(0) 
      if(verbose) cat("Warning: No outliers found in objNull.\n")
  } else {
      # FIX 2025-01-31: C++ kernel expects 1-based indices (it performs -1 internally).
      # objNull stores 0-based indices. We must convert to 1-based here.
      posOutlier <- as.integer(outlier_list$posOutlier + 1)
  }

  if (verbose) {
     cat("Samples in Null Model:", length(objNull_sample_ids), "\n")
     cat("Outliers:", length(posOutlier), "\n")
  }

  # 3. Match with File Samples
  # Call C++ helper to read IDs from gz file
  file_sample_ids <- read_ukb_sample_ids_cpp(dosage_file)
  
  match_result <- match(objNull_sample_ids, file_sample_ids)
  if (any(is.na(match_result))) {
      missing_count <- sum(is.na(match_result))
      if(missing_count == length(objNull_sample_ids)) {
          stop("NO samples matched between Null Model and Genotype file. Check ID formats.")
      }
      warning(paste(missing_count, "samples from Null Model not found in Genotype file. They will be ignored."))
  }
  
  # Logic: Filter `R` (residuals) to only include samples that exist in the file.
  valid_indices_in_objNull <- which(!is.na(match_result))
  R_subset <- R[valid_indices_in_objNull]
  
  # Indices in FILE (0-based) corresponding to the valid samples
  file_indices_0based <- as.integer(match_result[valid_indices_in_objNull] - 1)
  
  # Re-map outliers
  outlier_indices_1based = posOutlier + 1
  new_outlier_indices = match(outlier_indices_1based, valid_indices_in_objNull)
  # Filter out outliers that were dropped (not in file)
  new_outlier_indices = new_outlier_indices[!is.na(new_outlier_indices)]
  posOutlier_subset = as.integer(new_outlier_indices - 1)
  
  if (verbose) {
      cat("Matched Samples:", length(R_subset), "\n")
      cat("Mapped Outliers:", length(posOutlier_subset), "\n")
  }

  # 4. Load Phi Matrices
  # Helper to load and format phi matrix
  # Map IDs to indices in valid samples ONLY
  kept_ids = objNull_sample_ids[valid_indices_in_objNull]
  sample_map_idx <- setNames(seq_along(kept_ids), kept_ids)

  convert_phi <- function(phi_path) {
      # Read phi file
      p_data <- data.table::fread(phi_path, header = TRUE)
      
      if(nrow(p_data) == 0) return(matrix(0, 0, 3))
      
      # Assume columns: IID1, IID2, Phi_value (or similar names)
      colnames(p_data) <- c("i_id", "j_id", "phi_value")
      
      i_ids = as.character(p_data$i_id)
      j_ids = as.character(p_data$j_id)
      vals = as.numeric(p_data$phi_value)
      
      # Map to internal indices (1-based for R data.table)
      idx_i = sample_map_idx[i_ids]
      idx_j = sample_map_idx[j_ids]
      
      dt = data.table::data.table(i = idx_i, j = idx_j, val = vals)
      dt = dt[!is.na(i) & !is.na(j)] # Keep only pairs in our sample
      
      if(nrow(dt) == 0) return(matrix(0, 0, 3))
      
      # Return matrix (i, j, value)
      as.matrix(dt[, .(i, j, val)])
  }

  phi_A_mat <- convert_phi(paste0(phi_file_pattern, "A.txt"))
  phi_B_mat <- convert_phi(paste0(phi_file_pattern, "B.txt"))
  phi_C_mat <- convert_phi(paste0(phi_file_pattern, "C.txt"))
  phi_D_mat <- convert_phi(paste0(phi_file_pattern, "D.txt"))
  
  if (verbose) cat("Phi matrices loaded and mapped.\n")
  
  # 5. Get Total SNPs (Estimate)
  # Uses system call to support progress tracking
  total_snps = 0
  if (.Platform$OS.type == "unix") {
      count_cmd <- paste0("zcat \"", dosage_file, "\" | wc -l")
      total_snps = tryCatch({
          out = system(count_cmd, intern=TRUE, ignore.stderr=TRUE)
          as.integer(out) - 1 # header
      }, warning = function(w) 0, error = function(e) 0)
  }
  
  if(total_snps <= 0) total_snps = 10000000 # Dummy large number if unknown

  # 6. Run Analysis
  # Use default cutoff of 2.0 for SPA
  result <- SPAmixPlus_local_ukb_high_performance_streaming_cpp(
    geno_file = dosage_file,
    haplo_file = haplo_file, 
    output_file = output_file,
    file_match_idx = file_indices_0based,
    R_matched = R_subset,
    phi_A_mat = phi_A_mat,
    phi_B_mat = phi_B_mat,
    phi_C_mat = phi_C_mat,
    phi_D_mat = phi_D_mat,
    posOutlier = posOutlier_subset,
    total_snps = as.integer(total_snps),
    save_interval = as.integer(save_interval),
    cutoff = 2.0,  # Standard SPA cutoff
    MAF_cutoff = as.numeric(MAF_cutoff),
    MAC_cutoff = as.numeric(MAC_cutoff),
    verbose = as.logical(verbose)
  )
  
  return(result)
}
