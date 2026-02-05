#### test SPAmixPlus package using extdata from GRAB package
library(SPAmixPlus)

# --------------------------------------------------------------------------------

GRAB_extdata_dir = "/gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/"
sparse_GRM_file = paste0(GRAB_extdata_dir, "/SparseGRM.txt")
GenoFile = paste0(GRAB_extdata_dir, "/simuBGEN.bgen")
# GenoFileIndex = NULL # Auto-detect .bgi and .sample files
GenoFileIndex = c(paste0(GRAB_extdata_dir, "/simuBGEN.bgen.bgi"),
                  paste0(GRAB_extdata_dir, "/simuBGEN.sample"))
PhenoFile = paste0(GRAB_extdata_dir, "/simuPHENO.txt")

PhenoData = data.table::fread(PhenoFile)
head(PhenoData)


#### function to get null model residuals
SPA_G_Get_Resid = function(traits="survival/binary",
                           formula=NULL,
                           data=NULL,
                           pIDs=NULL,
                           gIDs=NULL,
                           range=c(-100,100),
                           length.out = 10000,
                           ...)
  
{
  if(traits=="survival"){
    Call = match.call()
    ### Fit a Cox model
    obj.coxph = coxph(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.coxph, range)
    
    ### Get the covariate matrix to adjust for genotype
    resid = obj.coxph$residuals
    
    re = resid
  }
  else if(traits=="binary"){
    Call = match.call()
    ### Fit a logistic model
    obj.logistic = glm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.logistic, range)
    
    ### Get the covariate matrix to adjust for genotype
    mu = obj.logistic$fitted.values
    resid = obj.logistic$y - mu
    re = resid
  }
  return(re)
}

check_input = function(pIDs, gIDs, obj, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = obj$residuals
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}

check_input_Resid = function(pIDs, gIDs, R, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = R
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}





#### Calculate p values ################################################################################################

resid = SPA_G_Get_Resid("binary",
                        BinaryPheno ~ AGE+GENDER+PC1+PC2,
                        data=PhenoData,
                        pIDs=PhenoData$IID,
                        gIDs=PhenoData$IID)



# step1: fit a null model

objNull_SPAmixPlus = SPAmixPlus.NullModel(resid ~ AGE+GENDER+PC1+PC2,
                                          data = PhenoData,
                                          subjData = IID, 
                                          sparseGRMFile = sparse_GRM_file,
                                          # method = "SPAmixPlus", 
                                          traitType = "Residual", 
                                          control = list(PC_columns = "PC1,PC2"))

#### calculate p-values from SPAmix+

# v1: do not estimate individual-specific in step0

output_dir = paste0("/gdata02/user/yuzhuoma/test_data/test_SPAmixPlus_package/test_GRAB_extdata/")
dir.create(output_dir, recursive = T)
OutputFile_SPAmixPlus_Global_v1 = paste0(output_dir,"/Output_SPAmixPlus_Global_v1_BGEN.txt")
file.remove(OutputFile_SPAmixPlus_Global_v1)
file.remove(paste0(OutputFile_SPAmixPlus_Global_v1, ".index"))

SPAmixPlus.Marker(objNull_SPAmixPlus,
                  GenoFile = GenoFile,
                  GenoFileIndex = GenoFileIndex,
                  OutputFile = OutputFile_SPAmixPlus_Global_v1,
                  control = list(AllMarkers = TRUE,
                                 min_maf_marker = 0.0001,
                                 outputColumns = "zScore",
                                 min_mac_marker = 1,
                                 AlleleOrder = "ref-first"))

Output_SPAmixPlus_Global_v1 = data.table::fread(OutputFile_SPAmixPlus_Global_v1)
head(Output_SPAmixPlus_Global_v1) # Note: The column BetaG denotes effect size estimate from SPAmix+

# v2: Estimate individual-specific in step0 to enhance computational efficieny in step 2
# =========================================================================
# Step 0: 预计算并建立 AF 估计 (只需运行一次)
# =========================================================================
# 无论以后要分析多少个表型，甚至分析不同的 SNP 子集，这步生成的 .bin 文件都是通用的。
# 它根据 GenoFile 里所有 SNP 的顺序，计算并存储了每个 SNP 的回归系数。

af_model_db = paste0(output_dir, "/SPAmixPlus_AF_Estimates_BGEN.bin")

# 控制参数：建议放宽过滤条件，让所有 SNP 都算出来存进去，Step 2 再去过滤。
control_step0 = list(min_maf_marker = 0.00001,  # 极低 MAF，保证几乎所有 SNP 都被记录
                     min_mac_marker = 1,
                     impute_method = "mean", 
                     omp_num_threads = 10,
                     AlleleOrder = "ref-first")      # Step 0 比较慢，可以用多线程加速

cat("启动 Step 0：构建全基因组 AF 模型库...\n")

SPAmixPlus.AF(GenoFile = GenoFile, 
              objNull = objNull_SPAmixPlus, 
              outputFile = af_model_db,        # 结果保存于此
              GenoFileIndex = GenoFileIndex,
              control = control_step0)

cat("Step 0 完成。数据库已保存。\n")


# inspect_step0_file(af_model_db)


# (可选) 检查文件内容确保不是全0
# inspect_step0_file(af_model_db, n_pcs = 4)

# =========================================================================
# Step 2: 关联分析 (针对特定表型)
# =========================================================================
# 这里的 control 是给 Step 2 筛选用的，可以比 Step 0 更严格。

OutputFile_SPAmixPlus_Global_v2 = paste0(output_dir,"/Output_SPAmixPlus_Global_v2_BGEN.txt")
file.remove(OutputFile_SPAmixPlus_Global_v2)
file.remove(paste0(OutputFile_SPAmixPlus_Global_v2, ".index"))

control_analysis = list(
  AllMarkers = TRUE,       # 分析所有标记
  min_maf_marker = 0.0001, # 可以在这里进行实际的分析过滤
  outputColumns = "zScore",
  afFile = af_model_db,     # <--- 关键点：传入二进制文件路径，开启极速模式
  AlleleOrder = "ref-first"
)

cat("启动 Step 2：极速关联分析...\n")

SPAmixPlus.Marker(objNull_SPAmixPlus,
                  GenoFile = GenoFile,
                  GenoFileIndex = GenoFileIndex,
                  OutputFile = OutputFile_SPAmixPlus_Global_v2,
                  control = control_analysis)

cat("分析结束。\n")

Output_SPAmixPlus_Global_v2 = data.table::fread(OutputFile_SPAmixPlus_Global_v2)
head(Output_SPAmixPlus_Global_v2) # Note: The column BetaG denotes effect size estimate from SPAmix+


# check difference between v1 and v2

sum(Output_SPAmixPlus_Global_v1$Pvalue - Output_SPAmixPlus_Global_v2$Pvalue)^2




# 查看 Step 0 二进制文件的工具函数
inspect_step0_file = function(bin_file, n_snps_to_read, n_pcs) {
  con = file(bin_file, "rb")
  
  # 每一个 SNP 的数据块大小 (以字节为单位)
  # 1 个 int (status) + (n_pcs + 1) 个 double (betas)
  record_size = 4 + (n_pcs + 1) * 8
  
  cat("Binary File Inspection for:", basename(bin_file), "\n")
  cat("Record Size:", record_size, "bytes per SNP\n\n")
  
  for(i in 1:n_snps_to_read) {
    if(isIncomplete(con)) break
    
    # 1. 读取状态 (Status)
    status = readBin(con, what = "integer", n = 1, size = 4)
    
    # 2. 读取系数 (Betas)
    betas = readBin(con, what = "numeric", n = n_pcs + 1, size = 8)
    
    # 翻译状态
    model_type = switch(as.character(status), 
                        "0" = "Mean (Fallback)", 
                        "1" = "Linear Model", 
                        "2" = "Logistic (Preferred)", 
                        "Unknown")
    
    cat(sprintf("SNP #%d: Type = %s\n", i, model_type))
    cat("  Betas (Intercept + PCs):", paste(round(betas, 4), collapse = ", "), "\n")
  }
  close(con)
}

# 使用示例 (请替换为你的实际文件路径)
# inspect_step0_file("/gdata02/user/yuzhuoma/test_data/.../SPAmixPlus_AF_Model_All.bin")

inspect_step0_file(af_model_db,
                   n_snps_to_read = nrow(Output_SPAmixPlus_Global_v2),
                   n_pcs = 2)


# check PLINK and BGEN difference

output_v2_PLINK = data.table::fread("/gdata02/user/yuzhuoma/test_data/test_SPAmixPlus_package/test_GRAB_extdata/Output_SPAmixPlus_Global_v2_PLINK.txt")
output_v2_BGEN = data.table::fread("/gdata02/user/yuzhuoma/test_data/test_SPAmixPlus_package/test_GRAB_extdata/Output_SPAmixPlus_Global_v2_BGEN.txt")


sum(output_v2_PLINK$Pvalue - output_v2_BGEN$Pvalue)^2
