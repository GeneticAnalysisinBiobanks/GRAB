# GRAB

### How to install and load this package

```{r}      
library(devtools)  # author version: 2.3.0
install_github("GeneticAnalysisinBiobanks/GRAB", INSTALL_opts=c("--no-multiarch"))  # The INSTALL_opts is required in Windows OS.
library(GRAB)
?GRAB.ReadGeno
?GRAB.SPACox
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
