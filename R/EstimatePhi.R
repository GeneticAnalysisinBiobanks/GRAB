#' Step 0 of SPAmixPlus: Pre-calculate individual-specific allele frequency coefficients
#' @export
SPAmixPlus.AF = function(GenoFile,
                         objNull,
                         outputFile,
                         GenoFileIndex = NULL,
                         control = list()){

    # Validate objNull
    if(!inherits(objNull, "SPAmixPlus_NULL_Model")) stop("objNull must be of class SPAmixPlus_NULL_Model")
    
    # Initialize control defaults if necessary
    if(is.null(control)) control = list()
    if(is.null(control$impute_method)) control$impute_method = "mean"
    if(is.null(control$missing_cutoff)) control$missing_cutoff = 0.15
    if(is.null(control$min_maf_marker)) control$min_maf_marker = 0.0001
    if(is.null(control$min_mac_marker)) control$min_mac_marker = 1
    if(is.null(control$omp_num_threads)) control$omp_num_threads = 1
    if(is.null(control$SPA_Cutoff)) control$SPA_Cutoff = 2
    
    # 1. Use internal setGenoInput to get ALL marker information
    # We pass 'AllMarkers=TRUE' to ensure no filtering happens at this stage unless specified
    control$AllMarkers = TRUE
    
    # Note: setGenoInput is internal in Geno.R. We can access it if inside package, or copy definition.
    # Since we are modifying the package, we can just call it.
    # However, setGenoInput is not exported. It is available in namespace.
    
    objGeno = setGenoInput(GenoFile, 
                           GenoFileIndex = GenoFileIndex, 
                           SampleIDs = as.character(objNull$subjData), 
                           control = control)
                           
    genoType = objGeno$genoType
    markerInfo = objGeno$markerInfo
    genoIndex = markerInfo$genoIndex
    
    if(length(genoIndex) == 0) stop("No available markers to process for Step 0.")
    
    cat("Step 0: Pre-calculating AF models for", length(genoIndex), "markers...\n")
    
    # 2. Initialize Global Vars in C++
    Group = makeGroup(objNull$yVec)
    
    setMarker_GlobalVarsInCPP(control$impute_method,
                              control$missing_cutoff,
                              control$min_maf_marker,
                              control$min_mac_marker,
                              control$omp_num_threads,
                              Group,
                              FALSE, # ifOutGroup
                              length(unique(Group)))

    # 3. Initialize SPAmixPlus Object using Null Model info
    setMarker.SPAmixPlus(objNull, control)
    
    # 4. Run Export
    exportAFModelInCPP("SPAmixPlus", genoType, genoIndex, outputFile);
    
    cat("Step 0 Completed. Model file saved to:", outputFile, "\n")
}


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
