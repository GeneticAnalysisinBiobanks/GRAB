


# step0
# 控制参数：建议放宽过滤条件，让所有 SNP 都算出来存进去，Step 2 再去过滤。
control_step0 = list(min_maf_marker = 0.00001,  # 极低 MAF，保证几乎所有 SNP 都被记录
                     min_mac_marker = 1,
                     impute_method = "mean", 
                     omp_num_threads = 10)      # Step 0 比较慢，可以用多线程加速

cat("启动 Step 0：构建全基因组 AF 模型库...\n")

SPAmixPlus.AF(GenoFile = GenoFile, 
              objNull = objNull_SPAmixPlus, 
              outputFile = af_model_db,        # 结果保存于此
              control = control_step0)

cat("Step 0 完成。数据库已保存。\n")

# step1
objNull_SPAmixPlus = SPAmixPlus.NullModel(resid ~ AGE+GENDER+PC1+PC2,
                                          data = PhenoData,
                                          subjData = IID, 
                                          sparseGRMFile = sparse_GRM_file,
                                          # method = "SPAmixPlus", 
                                          traitType = "Residual", 
                                          control = list(PC_columns = "PC1,PC2"))


# plink
OutputFile_SPAmixPlus_Global_v2 = paste0(output_dir,"/Output_SPAmixPlus_Global_v2_PLINK.txt")
file.remove(OutputFile_SPAmixPlus_Global_v2)
file.remove(paste0(OutputFile_SPAmixPlus_Global_v2, ".index"))

control_analysis = list(
  AllMarkers = TRUE,       # 分析所有标记
  min_maf_marker = 0.0001, # 可以在这里进行实际的分析过滤
  outputColumns = "zScore",
  afFile = af_model_db     # output from exportAFModelInCPP
)

cat("启动 Step 2：极速关联分析...\n")

SPAmixPlus.Marker(objNull_SPAmixPlus,
                  GenoFile = GenoFile,
                  OutputFile = OutputFile_SPAmixPlus_Global_v2,
                  control = control_analysis)





                