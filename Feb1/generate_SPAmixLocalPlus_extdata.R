output_dir = "/gdata02/user/yuzhuoma/GRAB_20260202/GRAB-main/inst/extdata/SPAmixPlus_extdata/"

library(data.table)
library(dplyr)
library(tidyr)

# Ancestry 1 Files
Anc1_Hap_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry1_HapCount.txt"
Anc1_Dos_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry1_Dosage.txt"

# Ancestry 2 Files
Anc2_Hap_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry2_HapCount.txt"
Anc2_Dos_File = "/gdata02/user/yuzhuoma/proj1/local_ance/data/global_PCs/genotype_rep_nSNP_50000/rep_1/PC_Data_Ancestry2_Dosage.txt"






Ancestry1_HapCount_temp = data.table::fread(Anc1_Hap_File)
Ancestry1_Dosage_temp = data.table::fread(Anc1_Dos_File)
Ancestry2_HapCount_temp = data.table::fread(Anc2_Hap_File)
Ancestry2_Dosage_temp = data.table::fread(Anc2_Dos_File)


Ancestry1_HapCount = Ancestry1_HapCount_temp[1:1000, 1:1001]
Ancestry1_Dosage = Ancestry1_Dosage_temp[1:1000, 1:1001]
Ancestry2_HapCount = Ancestry2_HapCount_temp[1:1000, 1:1001]
Ancestry2_Dosage = Ancestry2_Dosage_temp[1:1000, 1:1001]



setwd(output_dir)
data.table::fwrite(Ancestry1_HapCount, file = "Ancestry1_HapCount.txt")
data.table::fwrite(Ancestry2_HapCount, file = "Ancestry2_HapCount.txt")
data.table::fwrite(Ancestry1_Dosage, file = "Ancestry1_Dosage.txt")
data.table::fwrite(Ancestry2_Dosage, file = "Ancestry2_Dosage.txt")




