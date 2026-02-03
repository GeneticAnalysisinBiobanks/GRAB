

# ============================== SPAmixLocalPlus: Step0 =================================

#' Estimate Ancestry-Specific Kinship Coefficient (Phi)
#' 
#' @description
#' Estimates the ancestry-specific kinship coefficient (Phi) for admixed populations.
#' This function supports multiple modes of operation including in-memory matrix processing and 
#' streaming processing for large datasets like UK Biobank.
#' 
#' @param GRM Three-column GRM (data.frame, matrix, or file path). Columns: ID1, ID2, Value.
#' @param AncestryData A list where each element represents an ancestry (ancestry_id) and contains genetic data configuration.
#' @param SampleIDs Vector of sample IDs (required).
#' @param OutputDir Directory to save results.
#' @param Scenarios Vector of scenarios to calculate (default: c("A", "B", "C", "D")).
#' @param Threshold Phi threshold (default 0.0).
#' @param TaskID Task ID for parallel jobs (default 0). Used in output filename.
#' @param MAF_cutoff MAF cutoff for filtering (default 0.2, matching UKB protocol).
#' @param Mode "memory" (matrix inputs), "stream_ukb" (UKB standard files), "stream_file" (general files read into memory).
#' @param InputFormat For file inputs: "MemoryMatrix" (Rows=Samples, Cols=SNPs, read via fread) or "UKB" (Rows=SNPs, Cols=Samples).
#' @param ... Additional arguments.
#' 
#' @details
#' This function is crucial for account for population structure and relatedness in admixed populations.
#' 
#' ## Modes
#' \itemize{
#'   \item \code{memory}: Suitable for smaller datasets where genotype matrices can be loaded into memory.
#'   \item \code{stream_ukb}: Optimized for UK Biobank file structures.
#'   \item \code{stream_file}: Supports custom file paths for streaming.
#' }
#' 
#' @return No return value. Results are saved to \code{OutputDir}.
#' 
#' @examples
#' \dontrun{
#' # SPAmixPlus.EstimatePhi(GRM = grm_df, 
#' #                        AncestryData = list(1 = list(ukb_dir = "path")),
#' #                        SampleIDs = sample_ids,
#' #                        OutputDir = "output")
#' }
#' 
#' @import data.table
#' @export
SPAmixPlus.EstimatePhi = function(GRM, 
                                  AncestryData, 
                                  SampleIDs, 
                                  OutputDir, 
                                  Scenarios = c("A", "B", "C", "D"),
                                  Threshold = 0.0,
                                  TaskID = 0,
                                  MAF_cutoff = 0.2,
                                  Mode = NULL,
                                  InputFormat = "MemoryMatrix")
{
  # 1. Setup Output Directory
  if(!dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
  
  # 2. Process GRM (Inputs: File path or Object)
  if(is.character(GRM)) {
    if(!file.exists(GRM)) stop("GRM file not found.")
    grm_dt = data.table::fread(GRM)
  } else {
    grm_dt = data.table::as.data.table(GRM)
  }
  
  # Standardize GRM columns (ID1, ID2, Value)
  if(ncol(grm_dt) < 3) stop("GRM must have at least 3 columns (ID1, ID2, Value).")
  if(!all(c("ID1", "ID2") %in% names(grm_dt))) {
    colnames(grm_dt)[1:3] = c("ID1", "ID2", "Value")
  }
  
  # Map IDs to Indices (0-based for C++)
  sample_map = match(SampleIDs, SampleIDs) - 1 
  names(sample_map) = SampleIDs
  
  grm_dt$idx1 = match(grm_dt$ID1, SampleIDs) - 1
  grm_dt$idx2 = match(grm_dt$ID2, SampleIDs) - 1
  
  valid_grm = grm_dt[!is.na(idx1) & !is.na(idx2)]
  if(nrow(valid_grm) < nrow(grm_dt)) warning("Some GRM pairs contain IDs not in SampleIDs list.")
  
  # 3. Determine Execution Mode
  if(is.null(Mode)) {
     first_anc = AncestryData[[1]]
     # Heuristic: If 'ukb_dir' is present, use stream_ukb.
     # If 'hapcount' is a file path, check InputFormat.
     if(!is.null(first_anc$ukb_dir)) {
        Mode = "stream_ukb"
     } else if(is.character(first_anc$hapcount) || (is.list(first_anc) && is.character(first_anc[[1]]))) {
        # File paths provided directly
        Mode = "file_matrix" 
     } else {
        Mode = "memory"
     }
  }
  
  cat(paste0("Running EstimatePhi in [", Mode, "] mode.\n"))
  
  anc_names = names(AncestryData)
  if(is.null(anc_names)) anc_names = as.character(seq_along(AncestryData))
  
  for(i_anc in seq_along(AncestryData)) {
    anc_id = anc_names[i_anc]
    anc_dat = AncestryData[[i_anc]]
    
    for(scenario in Scenarios) {
      cat(paste0("  Processing Ancestry ", anc_id, " | Scenario ", scenario, " ...\n"))
      
      Phi_Values = NULL
      
      if(Mode == "stream_ukb") {
        # UKB Mode: Construct file lists and pass to C++ stream kernel
        if(is.null(anc_dat$ukb_dir) || is.null(anc_dat$bim_file)) {
          stop("For 'stream_ukb' mode, provide 'ukb_dir' and 'bim_file'.")
        }
        
        chroms = anc_dat$chromosomes
        if(is.null(chroms)) chroms = 1:22
        
        # Determine strict paths based on UKB naming convention
        # We construct the exact list of files to pass to C++
        hap_files = file.path(anc_dat$ukb_dir, paste0("ukb22418_c", chroms, ".anc", anc_id, ".hapcount.txt.gz"))
        dos_files = file.path(anc_dat$ukb_dir, paste0("ukb22418_c", chroms, ".anc", anc_id, ".dosage.txt.gz"))
        
        # Verify files exist? C++ will warn/skip if missing.
        
        res = estimatePhiStreamCPP(hap_files = hap_files,
                                   dos_files = dos_files,
                                   bim_file = anc_dat$bim_file,
                                   ancestry_id = anc_id,
                                   scenario = scenario,
                                   pair_idx1 = valid_grm$idx1,
                                   pair_idx2 = valid_grm$idx2,
                                   chromosomes = chroms,
                                   phi_threshold = Threshold,
                                   maf_cutoff = MAF_cutoff)
        
        Phi_Values = res$phi_values
        
      } else if (Mode == "file_matrix") {
        # General File Mode: Read file -> Memory Matrix -> Compute
        # Supports "Simulation" format (.txt/gz, Rows=Samples, Header)
        
        hap_file = anc_dat$hapcount
        dos_file = anc_dat$dosage
        
        cat("    Reading Dosage file... ")
        # InputFormat="MemoryMatrix" implies Rows=Samples. C++ needs Cols=Samples.
        # So we read and Transpose.
        
        # Assuming header check:
        # If user says InputFormat="MemoryMatrix", we treat it as standard simulation format
        
        dos_dt = data.table::fread(dos_file, header = TRUE)
        # Check alignment: First col is ID?
        # Simulation script: first col is ID.
        # We should verify IDs match SampleIDs or just assume order.
        # Safest: Use SampleIDs to subset/reorder if row names present.
        
        # Drop first column (ID) and convert to matrix
        ids_in_file = dos_dt[[1]]
        dos_mat = as.matrix(dos_dt[, -1])
        # Transpose to SNPs x Samples
        dos_mat = t(dos_mat)
        
        cat("Reading HapCount file...\n")
        hap_dt = data.table::fread(hap_file, header = TRUE)
        hap_mat = as.matrix(hap_dt[, -1])
        hap_mat = t(hap_mat)
        
        # Call Memory Kernel
        res = computePhiRatiosCPP(hapcount_matrix = hap_mat,
                                  dosage_matrix = dos_mat,
                                  pair_idx1 = valid_grm$idx1,
                                  pair_idx2 = valid_grm$idx2,
                                  scenario = scenario,
                                  phi_threshold = Threshold,
                                  maf_cutoff = MAF_cutoff)
        
        if(length(res$valid_counts) > 0) {
             Phi_Values = res$ratio_sums / ifelse(res$valid_counts==0, 1, res$valid_counts)
             Phi_Values[res$valid_counts == 0] = NA
        }
        
      } else {
        # Memory Mode (Lists provided directly)
        hap = anc_dat$hapcount
        dos = anc_dat$dosage
        
        # Align orientation if needed (Samples x SNPs -> SNPs x Samples)
        if(nrow(hap) == length(SampleIDs) && ncol(hap) != length(SampleIDs)) {
          hap = t(hap)
          dos = t(dos)
        }
        
        res = computePhiRatiosCPP(hapcount_matrix = hap,
                                  dosage_matrix = dos,
                                  pair_idx1 = valid_grm$idx1,
                                  pair_idx2 = valid_grm$idx2,
                                  scenario = scenario,
                                  phi_threshold = Threshold,
                                  maf_cutoff = MAF_cutoff)
        
        if(length(res$valid_counts) > 0) {
             Phi_Values = res$ratio_sums / ifelse(res$valid_counts==0, 1, res$valid_counts)
             Phi_Values[res$valid_counts == 0] = NA
        }
      }
      
      # ---------------------------------------------------------
      # SAVE OUTPUT
      # ---------------------------------------------------------
      if(!is.null(Phi_Values)) {
         # Ensure Phi_Values is a standard vector (handle 1-col matrix from C++)
         Phi_Values = as.vector(Phi_Values)

         is_diff = valid_grm$idx1 != valid_grm$idx2
         idx1_vec = valid_grm$idx1[is_diff]
         idx2_vec = valid_grm$idx2[is_diff]
         
         # Expand order: i->j, j->i
         final_i_idx = as.vector(rbind(idx1_vec, idx2_vec))
         final_j_idx = as.vector(rbind(idx2_vec, idx1_vec))
         
         final_ids_i = SampleIDs[final_i_idx + 1]
         final_ids_j = SampleIDs[final_j_idx + 1]
         
         res_df = data.table::data.table(
           i = final_ids_i,
           j = final_ids_j,
           phi_value = Phi_Values
         )
         
         res_df = res_df[!is.na(phi_value)]
         if(Threshold > 0) res_df = res_df[phi_value >= Threshold]
         
         out_filename = sprintf("phi_result_ancestry%s_scenario%s_task%s.txt", 
                                anc_id, scenario, sprintf("%03d", as.integer(TaskID)))
         out_path = file.path(OutputDir, out_filename)
         
         data.table::fwrite(res_df, out_path, sep="\t", quote=FALSE)
         cat(paste0("    Saved: ", out_filename, "\n"))
      }
    }
  }
  
  cat("Phi Estimation Completed.\n")
}



# ============================== SPAmixLocalPlus: Step2 =================================

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
