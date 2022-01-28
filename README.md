# GRAB

### How to install and load this package

Rtools (https://cran.rstudio.com/bin/windows/Rtools/) should be installed before installing this package.

```{r}      
library(devtools)  # author version: 2.3.0, use install.packages("devtools") first
install_github("GeneticAnalysisinBiobanks/GRAB", INSTALL_opts=c("--no-multiarch"))  # The INSTALL_opts is required in Windows OS.
library(GRAB)
```

### Replicate the genotype simulation in the package
```{r}   
## Commen Variants with MAF ranging from (0.05, 0.5)
set.seed(12345)
OutList = GRAB.SimuGMat(nSub = 500, nFam = 50, FamMode = "10-members", nSNP = 10000,
                        MaxMAF = 0.5, MinMAF = 0.05)
GenoMat = OutList$GenoMat 
MissingRate = 0.05
indexMissing = sample(length(GenoMat), MissingRate * length(GenoMat))
GenoMat[indexMissing] = -9
extDir = system.file("extdata", package = "GRAB")
extPrefix = paste0(extDir, "/simuPLINK")
GRAB.makePlink(GenoMat, extPrefix)

setwd(extDir)
system("plink --file simuPLINK --make-bed --out simuPLINK")
system("plink --bfile simuPLINK --recode A --out simuRAW")
system("plink2 --bfile simuPLINK --export bgen-1.2 bits=8 ref-first --out simuBGEN  # UK Biobank use 'ref-first'")
system("bgenix -g simuBGEN.bgen -index")

## Rare Variants with MAF ranging from (0.0001, 0.01)
set.seed(34567)
OutList = GRAB.SimuGMat(nSub = 500, nFam = 50, FamMode = "10-members", nSNP = 10000,
                        MaxMAF = 0.01, MinMAF = 0.0001)
GenoMat = OutList$GenoMat 
MissingRate = 0.05
indexMissing = sample(length(GenoMat), MissingRate * length(GenoMat))
GenoMat[indexMissing] = -9
extDir = system.file("extdata", package = "GRAB")
extPrefix = paste0(extDir, "/simuPLINK_RV")
GRAB.makePlink(GenoMat, extPrefix)

setwd(extDir)
system("plink --file simuPLINK_RV --make-bed --out simuPLINK_RV")
system("plink --bfile simuPLINK_RV --recode A --out simuRAW_RV")
system("plink2 --bfile simuPLINK_RV --export bgen-1.2 bits=8 ref-first --out simuBGEN_RV  # UK Biobank use 'ref-first'")
system("bgenix -g simuBGEN_RV.bgen -index")
```

### Replicate the sparse GRM generation in the package
```{r}
GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
PlinkFile = tools::file_path_sans_ext(GenoFile)  # remove file extension
nPartsGRM = 2
for(partParallel in 1:nPartsGRM){
  getTempFilesFullGRM(PlinkFile, nPartsGRM, partParallel) # this function only supports Linux
}
SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
getSparseGRM(PlinkFile, nPartsGRM, SparseGRMFile)
data.table::fread(SparseGRMFile)
#            ID1      ID2     Value
#    1:     f1_1     f1_1 0.9591102
#    2:     f1_2     f1_2 1.0227757
#    3:     f1_3     f1_3 1.0278691
#    4:     f1_4     f1_4 1.0138606
#    5:     f1_5     f1_1 0.4640743
#   ---                            
# 2546: Subj-496 Subj-496 1.0031059
# 2547: Subj-497 Subj-497 0.9952855
# 2548: Subj-498 Subj-498 0.9829549
# 2549: Subj-499 Subj-499 1.0086861
# 2550: Subj-500 Subj-500 0.9786987
```

### 


### Replicate the phenotype simulation in the package
```{r}
FamFile = system.file("extdata", "simuPLINK.fam", package = "GRAB")
FamData = read.table(FamFile)
IID = FamData$V2  # Individual ID
n = length(IID)
set.seed(678910)
## The below is just to demonstrate the functions of GRAB package
Pheno = data.frame(IID = IID, Cova1 = rnorm(n), Cova2 = rbinom(n, 1, 0.5), 
                   binary = rbinom(n, 1, 0.5),
                   ordinal = rbinom(n, 3, 0.3),
                   quantitative = rnorm(n),
                   time = runif(n),
                   event = rbinom(n, 1, 0.2))
PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
write.table(Pheno, PhenoFile, row.names = F, quote = F, sep = "\t")
```

### Replicate the regionFile in the package
```{r}
set.seed(101112)
RegionData = data.frame(REGION = paste0("Region_", rep(1:100,each=100)),
                        MARKER = paste0("SNP_",1:10000),
                        ANNO1 = rbinom(10000, 1, 0.5),
                        ANNO2 = runif(10000))
```

### Step 1: Fit a null model using the sparse GRM and phenotype data
```{r}
GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
# GenoFile = system.file("extdata", "simuBGEN.bgen", package = "GRAB")  # BGEN file input is also supported in null model fitting
PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
Pheno = read.table(PhenoFile, header = T)
SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
objNullFile = system.file("results", "objNull.RData", package = "GRAB")

objNull = GRAB.NullModel(as.factor(ordinal) ~ Cova1 + Cova2, 
                         data = Pheno, 
                         subset = (event==0), 
                         subjData = Pheno$IID, 
                         method = "POLMM", 
                         traitType = "ordinal", 
                         GenoFile = GenoFile,
                         SparseGRMFile = SparseGRMFile)

objNull$tau    # 0.5655751
                         
save(objNull, file = objNullFile)
```

### Step 2a: Perform a marker-level association analysis
```{r}
objNullFile = system.file("results", "objNull.RData", package = "GRAB")
load(objNullFile)

OutputDir = system.file("results", package = "GRAB")
OutputFile = paste0(OutputDir, "/simuMarkerOutput.txt")
GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")

# If 'OutputFile' and 'OutputFileIndex' have existed, users might need to remove them first.
GRAB.Marker(objNull, 
            GenoFile = GenoFile,
            OutputFile = OutputFile)
            
## Check 'OutputFile' for the analysis results
```

### Step 2b: Perform a genome-wide region-level association analysis
```{r}
objNullFile = system.file("results", "objNull.RData", package = "GRAB")
load(objNullFile)

OutputDir = system.file("results", package = "GRAB")
OutputFile = paste0(OutputDir, "/simuRegionOutput.txt")
GenoFile = system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
RegionFile = system.file("extdata", "simuRegion.txt", package = "GRAB")
SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")

## make sure the output files does not exist at first
file.remove(OutputFile)
file.remove(paste0(OutputFile, ".markerInfo"))
file.remove(paste0(OutputFile, ".index"))

GRAB.Region(objNull, 
            GenoFile = GenoFile,
            OutputFile = OutputFile,
            RegionFile = RegionFile,
            SparseGRMFile = SparseGRMFile)
            
data.table::fread(OutputFile)

## additional columns of "zScore", "nSamplesInGroup", "AltCountsInGroup", "AltFreqInGroup"
## We do not recommend adding too many columns for all markers

file.remove(OutputFile)
file.remove(paste0(OutputFile, ".markerInfo"))
file.remove(paste0(OutputFile, ".index"))
GRAB.Region(objNull, 
            GenoFile = GenoFile,
            OutputFile = OutputFile,
            RegionFile = RegionFile,
            RegionAnnoHeader = c("ANNO1", "ANNO2"),
            SparseGRMFile = SparseGRMFile,
            control = list(outputColumns = c("beta", "seBeta", "zScore","nSamplesInGroup","AltCountsInGroup","AltFreqInGroup")))
            
data.table::fread(OutputFile)
            
## Check 'OutputFile' for the analysis results
```

### NOTE for BGEN file input

The current version of GRAB package only supports BGEN v1.2 using 8 bits compression (faster than using 16 bits). For example, if plink2 is used to make BGEN file, please refer to https://www.cog-genomics.org/plink/2.0/data#export
```
plink2 --bfile input --out output --export bgen-1.2 bits=8
```

The index file for BGEN file is also required (filename extension is ".bgen.bgi"). Please refer to https://enkre.net/cgi-bin/code/bgen/wiki/bgenix 


### More detailed information for package developers

For SPACox method, we should create the following files

```{r}
/src/SPACox.cpp
/src/SPACox.hpp
/R/SPACox.R
```
and modify the following codes

```{r}
/src/Main.cpp  
# static SPACox::SPACoxClass* ptr_gSPACoxobj = NULL;
# void setSPACoxobjInCPP()
/R/GRAB_Null_Model.R
/R/GRAB_Marker.R
/R/GRAB_Region.R
```
