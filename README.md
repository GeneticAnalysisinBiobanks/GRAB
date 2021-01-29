# GRAB

### How to install and load this package

```{r}      
library(devtools)  # author version: 2.3.0
install_github("GeneticAnalysisinBiobanks/GRAB")
library(GRAB)
?GRAB::readGeno  
GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
GenoMat = GRAB.ReadGeno(GenoFile)
head(GenoMat)
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
