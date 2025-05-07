# GRAB

## Note on 2022-08-26

We make a website for GRAB package: https://wenjianbi.github.io/grab.github.io/. Please check it for more recent update. 

## Previous README before 2022-08-26

### How to install and load this package

Rtools (https://cran.rstudio.com/bin/windows/Rtools/) should be installed before installing this package.

```{r}      
library(remotes)  # remotes library requires less dependency packages than devtools
install_github("GeneticAnalysisinBiobanks/GRAB", INSTALL_opts=c("--no-multiarch"), ref="main")  # The INSTALL_opts is required in Windows OS.
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
data.table::fread(OutputFile)
#           Marker        Info   AltFreq AltCounts MissingRate     Pvalue        beta     seBeta
#     1:     SNP_1     1:1:G:A 0.3934316       587  0.04481434 0.77381083  0.03200175  0.1113516
#     2:     SNP_2     1:2:G:A 0.4561995       677  0.04993598 0.30789980 -0.10816640 -0.1060831
#     3:     SNP_3     1:3:G:A 0.3899329       581  0.04609475 0.77397169  0.03153778  0.1098174
#     4:     SNP_4     1:4:G:A 0.4433834       650  0.06145967 0.33178713  0.10265438  0.1057725
#     5:     SNP_5     1:5:G:A 0.2355705       351  0.04609475 0.59348431  0.06879171  0.1288732
#    ---                                                                                        
#  9996:  SNP_9996  1:9996:G:A 0.3424566       513  0.04097311 0.04678621 -0.21914094 -0.1102191
#  9997:  SNP_9997  1:9997:G:A 0.2686062       397  0.05377721 0.84223235 -0.02476919 -0.1244440
#  9998:  SNP_9998  1:9998:G:A 0.0669344       100  0.04353393 0.75590206 -0.06695376 -0.2153778
#  9999:  SNP_9999  1:9999:G:A 0.2377384       349  0.06017926 0.45544541 -0.09688906 -0.1298141
# 10000: SNP_10000 1:10000:G:A 0.1361186       202  0.04993598 0.37979954 -0.13401746 -0.1525933

# The below is to show some parameters in control list. For more details, please check the Details section in ?GRAB.Marker
GRAB.Marker(objNull, 
            GenoFile = GenoFile,
            OutputFile = OutputFile,
            control = list(nMarkersEachChunk = 1000,
                           ImputeMethod = "mean",
                           MissingRateCutoff = 0.05,
                           MinMAFCutoff = 0.1))

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
