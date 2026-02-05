#' SPAmixLocalPlus: Local Analysis for Related Samples in UK Biobank
#' 
#' @description
#' Implements the SPAmixPlus local analysis algorithm for related samples in UK Biobank data.
#' This function integrates the "local" analysis approach for high-performance streaming of BGEN/Dosage data.
#' It is designed to be a direct replacement for the `SPAmixPlus_local_related_batch_rcpp` function 
#' found in the standalone scripts.
#' 
#' @param geno_file Character. Path to UKB genotype dosage file (.txt.gz).
#' @param haplo_file Character. Path to UKB haplotype hapcount file (.txt.gz).
#' @param phi_dir Character. Path to directory containing phi (kinship) data.
#' @param objNull Object of class `SPAmixPlus_NULL_Model` returned by `SPAmixPlus.NullModel`.
#' @param ancestry_name Character. Ancestry ID (e.g., "1", "2").
#' @param pheno_idx Integer. Phenotype index to analyze (default 1).
#' @param save_interval Integer. Number of SNPs processed before saving/checkpointing (default 100).
#' @param cutoff Numeric. SPA cutoff threshold (default 2.0).
#' @param MAF_cutoff Numeric. Minor Allele Frequency cutoff (default 0.00001).
#' @param MAC_cutoff Numeric. Minor Allele Count cutoff (default 1).
#' @param output_file Character. Path to output result file.
#' @param verbose Logical. If TRUE, output detailed progress (default TRUE).
#' @details 
#' The \code{SPAmixLocalPlus} function helps users perform local ancestry analysis efficiently. 
#' It reads genotype dosage and haplotype count files in a streaming fashion, minimizing memory usage.
#' 
#' ## Input Files
#' \itemize{
#'   \item \code{geno_file}: Dosage file (likely .gz compressed).
#'   \item \code{haplo_file}: Haplotype count file (likely .gz compressed).
#' }
#' 
#' ## Phi Directory
#' The \code{phi_dir} must contain the pre-calculated Phi (kinship) matrices or files.
#' These are typically generated using `SPAmixPlus.EstimatePhi`.
#' 
#' ## Output
#' The function writes results directly to \code{output_file}. 
#' The output format includes columns for Marker ID, P-value, and other statistics.
#' 
#' @return No return value. Results are saved to \code{output_file}.
#' 
#' @examples
#' \dontrun{
#' # Example usage:
#' SPAmixLocalPlus(
#'   geno_file = "path/to/geno.txt.gz",
#'   haplo_file = "path/to/haplo.txt.gz",
#'   phi_dir = "path/to/phi_dir",
#'   objNull = objNull,
#'   ancestry_name = "1",
#'   output_file = "output.txt"
#' )
#' }
#' 
#' @export
SPAmixLocalPlus = function(geno_file,
                           haplo_file, 
                           phi_dir,
                           objNull,
                           ancestry_name,
                           pheno_idx = 1,
                           save_interval = 100,
                           cutoff = 2.0,
                           MAF_cutoff = 0.00001,
                           MAC_cutoff = 1,
                           output_file,
                           verbose = TRUE) {
  
  # 1. Validation
  if (verbose) cat("=== SPAmixLocalPlus (Integrated) ===\n")
  if (!file.exists(geno_file)) stop("Genotype file not found: ", geno_file)
  if (!file.exists(haplo_file)) stop("Haplotype file not found: ", haplo_file)
  if (!dir.exists(phi_dir)) stop("Phi directory not found: ", phi_dir)
  
  # Check class of objNull
  # Note: The user might call it standard list or class. 
  # We check for basic required components if class check fails or behaves weirdly.
  if (!inherits(objNull, "SPAmixPlus_NULL_Model") && !inherits(objNull, "list")) {
      warning("objNull is not of class 'SPAmixPlus_NULL_Model' or 'list'. Proceeding with caution.")
  }

  if (verbose) {
    cat("Genotype File:", geno_file, "\n")
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

  # 3. Match with UKB Samples
  # Call C++ helper to read IDs from gz file (using the passed geno_file)
  ukb_sample_ids <- read_ukb_sample_ids_cpp(geno_file)
  
  match_result <- match(objNull_sample_ids, ukb_sample_ids)
  if (any(is.na(match_result))) {
      missing_count <- sum(is.na(match_result))
      if(missing_count == length(objNull_sample_ids)) {
          stop("NO samples matched between Null Model and UKB Genotype file. Check ID formats.")
      }
      warning(paste(missing_count, "samples from Null Model not found in UKB Genotype file. They will be ignored."))
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

  # 4. Load Phi Matrix
  # Helper to load and format phi matrix
  # Map IDs to indices in valid samples ONLY
  kept_ids = objNull_sample_ids[valid_indices_in_objNull]
  sample_map_idx <- setNames(seq_along(kept_ids), kept_ids)

  convert_phi <- function(p_data) {
      if(length(p_data$phi_value) == 0) return(matrix(0, 0, 3))
      
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

  phi_A_data <- load_phi_data_ukb_cpp(phi_dir, ancestry_name, "A")
  phi_B_data <- load_phi_data_ukb_cpp(phi_dir, ancestry_name, "B")
  phi_C_data <- load_phi_data_ukb_cpp(phi_dir, ancestry_name, "C")
  phi_D_data <- load_phi_data_ukb_cpp(phi_dir, ancestry_name, "D")

  phi_A_mat <- convert_phi(phi_A_data)
  phi_B_mat <- convert_phi(phi_B_data)
  phi_C_mat <- convert_phi(phi_C_data)
  phi_D_mat <- convert_phi(phi_D_data)
  
  if (verbose) cat("Phi matrices loaded and mapped.\n")
  
  # 5. Get Total SNPs (Estimate)
  # Uses system call as in original script to support progress bar
  total_snps = 0
  if (.Platform$OS.type == "unix") {
      count_cmd <- paste0("zcat \"", geno_file, "\" | wc -l")
      total_snps = tryCatch({
          out = system(count_cmd, intern=TRUE, ignore.stderr=TRUE)
          as.integer(out) - 1 # header
      }, warning = function(w) 0, error = function(e) 0)
  }
  
  if(total_snps <= 0) total_snps = 10000000 # Dummy large number if unknown

  # 6. Run Analysis
  result <- SPAmixPlus_local_ukb_high_performance_streaming_cpp(
    geno_file = geno_file,
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
    cutoff = as.numeric(cutoff),
    MAF_cutoff = as.numeric(MAF_cutoff),
    MAC_cutoff = as.numeric(MAC_cutoff),
    verbose = as.logical(verbose)
  )
  
  return(result)
}
