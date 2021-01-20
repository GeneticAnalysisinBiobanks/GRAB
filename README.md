# GRAB

### How to install and load this package

```{r}      
library(devtools)  # author version: 2.3.0
install_github("GeneticAnalysisinBiobanks/GRAB")
library(GRAB)
?GRAB::readGeno  
GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
GenoMat = readGeno(GenoFile)
head(GenoMat)
```
