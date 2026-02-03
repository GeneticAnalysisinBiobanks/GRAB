
# step0

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


# step1
objNull_SPAmixPlus = SPAmixPlus.NullModel(resid ~ AGE+GENDER+PC1+PC2,
                                          data = PhenoData,
                                          subjData = IID, 
                                          sparseGRMFile = sparse_GRM_file,
                                          # method = "SPAmixPlus", 
                                          traitType = "Residual", 
                                          control = list(PC_columns = "PC1,PC2"))


# step2
  SPAmixLocalPlus(
    geno_file = geno_file_1,
    haplo_file = haplo_file_1,
    phi_dir = data_dir,          # Phi文件所在的目录
    objNull = objNull_SPAmixPlus,# Null Model 对象
    ancestry_name = "1",         # 对应文件名中的 Ancestry1
    pheno_idx = 1,               # 通常为 1
    save_interval = 100,         # 每计算多少个 SNP 保存一次
    MAF_cutoff = 0.00001,        # 过滤低频变异
    MAC_cutoff = 1,
    output_file = output_file_ance1,
    verbose = TRUE               # 显示进度
  )

