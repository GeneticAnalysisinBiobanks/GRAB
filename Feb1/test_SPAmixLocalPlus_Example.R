# ============================================================================
# Script: SPAmixLocalPlus_Example.R
# Purpose: Demonstrate how to use SPAmixLocalPlus for local ancestry analysis
# Data Source: /gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/SPAmixPlus_extdata/
# ============================================================================

# 1. 设置环境与加载包
library(data.table)
library(dplyr)
library(survival)
library(SPAmixPlus)

# 清理内存
gc()

# 2. 定义输入输出路径
# 数据所在目录 (根据您的要求)
data_dir <- "/gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/SPAmixPlus_extdata/"
GRAB_extdata_dir = "/gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/"

# 输出结果目录
output_dir <- "/gdata02/user/yuzhuoma/test_data/test_SPAmixPlus_package/test_local_example_output/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Phi (Ancestry Information) 文件所在目录 (这里直接指向 data_dir，因为 ls 显示 phi 结果也在里面)
# 注意：SPAmixLocalPlus 会自动在 phi_dir 下查找对应的 phi 文件
phi_dir <- data_dir  

# 3. 模拟或加载表型数据 (Phenotype Data)
# 由于您提供的 ls 列表中没有明确的 phenotype 文件，这里我们构建一个模拟的表型数据
# 假设样本ID与 Ancestry 文件中的一致 (f1_1 到 f5000_4 这种格式，或其他格式，需与基因型文件匹配)
# 我们先读取基因型文件的列名来获取样本 ID

# 读取一个基因型文件的头部来获取样本ID
temp_geno <- data.table::fread(paste0(data_dir, "Ancestry1_Dosage.txt"), nrows=0) 
# 注意：原文件是 Dosage.txt (Ind x SNP)，第一列是 IID
# 如果您使用的是转换后的 UKB 格式 (SNP x Ind)，则列名是样本ID (除了前5列)
# 这里为了稳健，假设我们需要构建 1000 个样本 (根据之前的上下文 nSNP=1000, nFam=5000等，这里我们实际读取一下)
# 为了演示，我们构造一个通用的流程

cat("正在读取样本列表...\n")
# 这里我们直接使用转换前的 Ind x SNP 文件获取 ID，或者如果您确信 ID 格式
All_IIDs <- data.table::fread(paste0(data_dir, "Ancestry1_Dosage.txt"), select = 1, header = T)$IID
n_samples <- length(All_IIDs)
cat(paste("检测到样本数:", n_samples, "\n"))


sparse_GRM_file = paste0(GRAB_extdata_dir, "/SparseGRM.txt")
PhenoFile = paste0(GRAB_extdata_dir, "/simuPHENO.txt")

PhenoData = data.table::fread(PhenoFile)
head(PhenoData)

PhenoData = PhenoData %>% mutate(IID = All_IIDs)
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



# 4. 拟合 Null Model (使用 SPAmixPlus.NullModel)
# 这是运行 Local Ancestry Test 的前置步骤
cat("正在拟合 Null Model (SPAmixPlus.NullModel)...\n")
objNull_SPAmixPlus = SPAmixPlus.NullModel(resid ~ AGE+GENDER+PC1+PC2,
                                          data = PhenoData,
                                          subjData = IID, 
                                          sparseGRMFile = sparse_GRM_file,
                                          # method = "SPAmixPlus", 
                                          traitType = "Residual", 
                                          control = list(PC_columns = "PC1,PC2"))






# 5. 定义输出文件路径
output_file_ance1 <- file.path(output_dir, "Result_Ancestry1.txt")
output_file_ance2 <- file.path(output_dir, "Result_Ancestry2.txt")

# 6. 运行 Local Ancestry Test - Ancestry 1
cat("\n=== 运行 Ancestry 1 分析 ===\n")
# 输入文件必须是之前转换好的 .txt.gz 格式 (UKB 格式: SNP x Ind)
geno_file_1  <- paste0(data_dir, "Ancestry1_Dosage_UKB.txt.gz")
haplo_file_1 <- paste0(data_dir, "Ancestry1_HapCount_UKB.txt.gz")

if(file.exists(geno_file_1) & file.exists(haplo_file_1)) {
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
  
  # 查看结果
  if(file.exists(output_file_ance1)){
    res1 <- data.table::fread(output_file_ance1)
    print(head(res1))
    cat("Ancestry 1 分析完成。显著结果 (P < 0.05):\n")
    print(res1 %>% filter(Pvalue < 0.05) %>% head())
  }
} else {
  cat("错误：Ancestry 1 的输入文件不存在，请检查之前的转换步骤。\n")
}

# 7. 运行 Local Ancestry Test - Ancestry 2
cat("\n=== 运行 Ancestry 2 分析 ===\n")
geno_file_2  <- paste0(data_dir, "Ancestry2_Dosage_UKB.txt.gz")
haplo_file_2 <- paste0(data_dir, "Ancestry2_HapCount_UKB.txt.gz")

if(file.exists(geno_file_2) & file.exists(haplo_file_2)) {
  SPAmixLocalPlus(
    geno_file = geno_file_2,
    haplo_file = haplo_file_2,
    phi_dir = data_dir,          
    objNull = objNull_SPAmixPlus,
    ancestry_name = "2",         # 对应文件名中的 Ancestry2
    pheno_idx = 1,
    save_interval = 100,
    MAF_cutoff = 0.00001,
    MAC_cutoff = 1,
    output_file = output_file_ance2,
    verbose = TRUE
  )
  
  # 查看结果
  if(file.exists(output_file_ance2)){
    res2 <- data.table::fread(output_file_ance2)
    print(head(res2))
  }
} else {
  cat("错误：Ancestry 2 的输入文件不存在。\n")
}

cat("\n所有任务已完成！结果保存于:", output_dir, "\n")