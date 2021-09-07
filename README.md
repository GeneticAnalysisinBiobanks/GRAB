# GRAB

### How to install and load this package

Rtools (https://cran.rstudio.com/bin/windows/Rtools/) should be installed before installing this package.

```{r}      
library(devtools)  # author version: 2.3.0, use install.packages("devtools") first
install_github("GeneticAnalysisinBiobanks/GRAB", INSTALL_opts=c("--no-multiarch"))  # The INSTALL_opts is required in Windows OS.
library(GRAB)
```

### For BGEN file input

The current version of GRAB package only supports BGEN v1.2 using 8 bits compression (faster than using 16 bits). For example, if plink2 is used to make BGEN file, please refer to https://www.cog-genomics.org/plink/2.0/data#export
```
plink2 --bfile input --out output --export bgen-1.2 bits=8
```

The index file for BGEN file is also required (filename extension is ".bgen.bgi"). Please refer to https://enkre.net/cgi-bin/code/bgen/wiki/bgenix 

```
theFolder/bgenix -g theFile.bgen -index
```

### Replicate the data simulation in the package
```{r}   
set.seed(12345)
OutList = GRAB.SimuGMatCommon(nSub = 500, nFam = 50, FamMode = "10-members", nSNP = 10000)
GenoMat = OutList$GenoMat 
MissingRate = 0.05
indexMissing = sample(length(GenoMat), MissingRate * length(GenoMat))
GenoMat[indexMissing] = -9
extDir = system.file("extdata", package = "GRAB")
extPrefix = paste0(extDir, "/simuPLINK")
makePlink(GenoMat, extPrefix)

setwd(extDir)
system("plink --file simuPLINK --make-bed --out simuPLINK")
system("plink --bfile simuPLINK --recode A --out simuRAW")
system("plink2 --bfile simuPLINK --export bgen-1.2 bits=8 ref-first --out simuBGEN  # UK Biobank use 'ref-first'")
system(bgenix -g simuBGEN.bgen --index)
```

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
