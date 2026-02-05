# ============================================================================
# 脚本名称: convert_to_txt.gz.R
# 功能: 将生成的 Ancestry*.txt (Individual x SNP) 转换为 UKB 格式 (SNP x Individual)
# 输出: .txt.gz 压缩文件，前5列为 CHROM, POS, ID, REF, ALT
# ============================================================================

library(data.table)
library(R.utils)

# 1. 设置工作目录 (您文件所在的目录)
work_dir <- "/gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/SPAmixPlus_extdata/"
if (!dir.exists(work_dir)) {
  stop(paste("目录不存在:", work_dir))
}
setwd(work_dir)
cat("当前工作目录:", getwd(), "\n")

# 2. 定义转换函数
convert_to_ukb_gz <- function(input_filename, output_filename) {
  
  start_time <- Sys.time()
  cat(paste0("\n=== 开始处理: ", input_filename, " ===\n"))
  
  if (!file.exists(input_filename)) {
    cat("错误: 文件不存在 ->", input_filename, "\n")
    return(FALSE)
  }
  
  # 读取数据 (假设第一列是 IID，后续列是 SNP)
  cat("  正在读取数据...\n")
  dt <- data.table::fread(input_filename, header = TRUE)
  
  # 获取样本 ID (第一列)
  sample_ids <- as.character(dt[[1]]) 
  n_samples <- length(sample_ids)
  
  # 获取 SNP ID (列名，排除第一列)
  snp_ids <- names(dt)[-1]
  n_snps <- length(snp_ids)
  
  cat(sprintf("  检测到数据维度: %d 个体 x %d SNPs\n", n_samples, n_snps))
  
  # 准备 SNP 信息列 (模拟信息)
  # 假设都在染色体 1，位置依次增加，REF=A, ALT=G
  ukb_df <- data.table(
    CHROM = rep(1, n_snps),
    POS   = 1000 + (1:n_snps) * 100,  # 模拟位置
    ID    = snp_ids,                  # 原始列名作为 SNP ID
    REF   = rep("A", n_snps),
    ALT   = rep("G", n_snps)
  )
  
  # 转置基因型矩阵 (从 Ind x SNP 转为 SNP x Ind)
  cat("  正在转置矩阵 (Ind x SNP -> SNP x Ind)...\n")
  geno_mat <- as.matrix(dt[, -1, with = FALSE]) # 去掉第一列 IID
  geno_mat_t <- t(geno_mat)                     # 转置
  
  # 将转置后的矩阵合并到 data.table
  # 注意：data.table 按列添加比较快
  cat("  合并为 UKB 格式表...\n")
  ukb_geno_dt <- as.data.table(geno_mat_t)
  names(ukb_geno_dt) <- sample_ids  # 列名设为样本ID
  
  # 最终合并：前5列 + 基因型列
  final_dt <- cbind(ukb_df, ukb_geno_dt)
  
  # 验证列数
  expected_cols <- 5 + n_samples
  if (ncol(final_dt) != expected_cols) {
    stop("列数计算错误，请检查代码逻辑。")
  }
  
  # 输出文件路径
  temp_txt <- gsub("\\.gz$", "", output_filename) # 如果输出名带.gz，先去掉
  if (!grepl("\\.txt$", temp_txt)) temp_txt <- paste0(temp_txt, ".txt") # 确保是 .txt
  
  # 写入 TXT
  cat(paste("  写入临时文件:", temp_txt, "...\n"))
  data.table::fwrite(final_dt, file = temp_txt, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 压缩为 GZ
  final_gz <- output_filename
  if (!grepl("\\.gz$", final_gz)) final_gz <- paste0(final_gz, ".gz")
  
  cat(paste("  压缩为:", final_gz, "...\n"))
  if (file.exists(final_gz)) file.remove(final_gz)
  R.utils::gzip(temp_txt, destname = final_gz, remove = TRUE) # remove=TRUE 删除源 .txt
  
  end_time <- Sys.time()
  cat(sprintf("  ✓ 完成。耗时: %.2f 秒\n", as.numeric(end_time - start_time, units = "secs")))
  return(TRUE)
}

# 3. 定义任务列表并执行
tasks <- list(
  list(input = "Ancestry1_HapCount.txt", output = "Ancestry1_HapCount_UKB.txt.gz"),
  list(input = "Ancestry2_HapCount.txt", output = "Ancestry2_HapCount_UKB.txt.gz"),
  list(input = "Ancestry1_Dosage.txt",   output = "Ancestry1_Dosage_UKB.txt.gz"),
  list(input = "Ancestry2_Dosage.txt",   output = "Ancestry2_Dosage_UKB.txt.gz")
)

cat("=== 开始批量转换任务 ===\n")

for (task in tasks) {
  tryCatch({
    convert_to_ukb_gz(task$input, task$output)
  }, error = function(e) {
    cat(paste("  !!! 任务失败:", task$input, "\n  错误信息:", e$message, "\n"))
  })
}

cat("\n所有任务执行完毕。\n")

# 4. 简单验证生成结果
cat("\n=== 验证生成文件 ===\n")
out_files <- sapply(tasks, function(x) x$output)
for (f in out_files) {
  if (file.exists(f)) {
    size_mb <- file.size(f) / 1024 / 1024
    cat(sprintf("[存在] %s (大小: %.2f MB)\n", f, size_mb))
    
    # 读取前几行验证列头
    tmp <- data.table::fread(f, nrows = 2)
    cols <- colnames(tmp)[1:5]
    if (all(cols == c("CHROM", "POS", "ID", "REF", "ALT"))) {
      cat("       -> 格式头验证通过 (CHROM, POS, ID, REF, ALT)\n")
    } else {
      cat("       -> [警告] 格式头不匹配！\n")
    }
  } else {
    cat(sprintf("[缺失] %s\n", f))
  }
}