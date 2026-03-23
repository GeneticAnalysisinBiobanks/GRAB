objNull <- readRDS(".test/hwe/maxFam1_TC_chr1.rds")
GRAB.Marker6(objNull, ".test/hwe/extract10.bed", ".test/hwe/marker6_extract10.tsv", overwrite = TRUE)
current_result <- fread(".test/hwe/marker6_extract10.tsv")
plink2_acount <- fread(".test/hwe/extract10.acount")
plink2_hardy <- fread(".test/hwe/extract10.hardy")

current_result[, c('Marker', 'AltCount', 'HWEpval')]
plink2_acount
plink2_hardy[, c('ID', 'HOM_A1_CT', 'HET_A1_CT', 'TWO_AX_CT', 'P')]
