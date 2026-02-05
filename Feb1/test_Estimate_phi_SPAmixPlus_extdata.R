
# =================================================================================================
# SPAmixPlus Phi Estimation Script (File-Based Version)
# =================================================================================================
# This script estimates Ancestry-Specific Phi using direct file inputs (.txt) without loading .RData.
# It uses the SPAmixPlus package function 'SPAmixPlus.EstimatePhi' in 'file_matrix' mode.
#
# SLURM Submission Examples:
# cd /path/to/script/
# sbatch --partition=bi2 -q low  -J phi_est --mem=100G -t 1-0:0 --array=1-3 -o log/%A_%a.log --wrap='Rscript test_SPAmixPlus_EstimatePhi_Migration.R $SLURM_ARRAY_TASK_ID'
# =================================================================================================


n.cpu = 1

# 3. Setup Logic for Thresholds (Based on your array task ID)
thresholdVec = c(0, 0.05, 0.1)
threshold = thresholdVec[n.cpu]
cat(paste0("Running Task ID: ", n.cpu, " | Threshold: ", threshold, "\n"))


# =================================================================================================
# 4. Input Configuration
# =================================================================================================

# 4.1 GRM Input
# -------------------------------------------------------------------------------------------------
# Path to the 3-column sparse GRM (ID1, ID2, Value).
# You can provide the file path string directly to the function.
sparse_GRM_file = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_nSNP_1e+05/sparse_GRM_3col.txt"

# 4.2 Sample ID List
# -------------------------------------------------------------------------------------------------
# Construct sample IDs (Must match the IDs or order in the input files).
nFam = 5000
IID_Vec = c()
for (FamID in c(1:nFam)) {
  IID_Vec_temp = paste0("f", FamID,"_", c(1:4))
  IID_Vec = c(IID_Vec, IID_Vec_temp)
}
cat(paste0("Total Samples: ", length(IID_Vec), "\n"))


# 4.3 Genetic Data Files (Direct Paths)
# -------------------------------------------------------------------------------------------------
# Define the paths to your text files.
# FORMAT ASSUMPTION: 
#   - Text files (.txt or .txt.gz)
#   - Rows represent SAMPLES, Columns represent SNPs (Standard simulation format)
#   - First column usually contains Sample IDs (The function will drop the first column automatically)
#   - Header exists

# Ancestry 1 Files
Anc1_Hap_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry1_HapCount.txt"
Anc1_Dos_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry1_Dosage.txt"

# Ancestry 2 Files
Anc2_Hap_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry2_HapCount.txt"
Anc2_Dos_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry2_Dosage.txt"

# 4.4 Ancestry Data Configuration List
# -------------------------------------------------------------------------------------------------
# Create a named list where keys ("1", "2") are Ancestry IDs.
# Each element is a list containing file paths for 'hapcount' and 'dosage'.
Ancestry_Config = list(
  "1" = list(
    hapcount = Anc1_Hap_File, 
    dosage = Anc1_Dos_File
  ),
  "2" = list(
    hapcount = Anc2_Hap_File, 
    dosage = Anc2_Dos_File
  )
)

# =================================================================================================
# 5. Run Estimation (SPAmixPlus.EstimatePhi)
# =================================================================================================

# Output Directory
# output_dir = "/gdata02/user/yuzhuoma/test_data/test_local/phi_results_v2_package"
# output_dir = "/gdata02/user/yuzhuoma/test_data/test_SPAmixPlus_package/local_UKB_test/phi_results_v30_package"

output_dir = "/gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/SPAmixPlus_extdata/"

dir.create(output_dir, recursive = T)

cat("Starting Phi Estimation...\n")

SPAmixPlus.EstimatePhi(
  GRM = sparse_GRM_file,            # Path to GRM file (Function will read it automatically)
  AncestryData = Ancestry_Config,   # The list of file paths created above
  SampleIDs = IID_Vec,              # Vector of Sample IDs for alignment and checking
  OutputDir = output_dir,           # Directory to save the phi results
  Scenarios = c("A", "B", "C", "D"),# The scenarios to calculate
  Threshold = threshold,            # Phi value cutoff (e.g., 0.05)
  TaskID = n.cpu,                   # Appended to filename (e.g., ..._task001.txt) to avoid overwrites
  
  # IMPORTANT MODE SETTING:
  Mode = "file_matrix",             # "file_matrix" tells the function to inputs are FILES
  # and that they need to be read (Samples x SNPs) and transposed
  # to match the internal computation format (SNPs x Samples).
  
  InputFormat = "MemoryMatrix"      # Confirms the file format is standard Simulation format 
  # (Rows=Samples, Header=True)
)

cat("Success! Phi estimation completed.\n")
cat("Results saved in:", output_dir, "\n")














