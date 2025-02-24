
# This function follows function getSparseGRM() 
getPairwiseIBD = function(PlinkFile,            # input path to PLINK file (without file extensions of bed/bim/fam).
                          SparseGRMFile,        # input path to SparseGRMFile from getSparseGRM() function.
                          PairwiseIBDFile,      # output path to save pairwise IBD to PairwiseIBDFile.
                          frqFile = NULL,       # input path to frq file corresponding to Plink file, default is PlinkFile.
                          tempDir = NULL,       # output path to save the temp files.
                          maxSampleNums = 2500, # read in at most 2500 subjects' genotypes for analysis.
                          minMafIBD = 0.01,     # Minimal value of MAF cutoff to select markers (default=0.01).
                          rm.tempFile = FALSE)  # a logical value indicating if the temp files will be deleted. (default=FALSE)
{
  cat("Noting that PlinkFile name here has no suffix (e.g. .bed or .frq).\n")
  
  bedFile = paste0(PlinkFile, ".bed")
  bimFile = paste0(PlinkFile, ".bim")
  famFile = paste0(PlinkFile, ".fam")
  
  if(is.null(frqFile)){
    frqFile = paste0(PlinkFile, ".frq")
  }else if(substr(frqFile, nchar(frqFile)-3, nchar(frqFile)) != ".frq"){
    frqFile = paste0(frqFile, ".frq")
  }
  if(!file.exists(frqFile))
    stop("Please check frqFile name ", frqFile, " is correct or this file exists!")
  
  if(is.null(tempDir))
    tempDir = system.file("PairwiseIBD", "temp", package = "GRAB")
  
  # read all genotype and pass to QC.
  GenoInfoMat = data.table::fread(frqFile)
  
  # read in the Sparse GRM.
  SparseGRMData = data.table::fread(SparseGRMFile)
  SparseGRMData$ID1 = as.character(SparseGRMData$ID1);
  SparseGRMData$ID2 = as.character(SparseGRMData$ID2)
  
  if(any(colnames(SparseGRMData) != c("ID1", "ID2", "Value")))
    stop("The column names of SparseGRMFile should be ['ID1', 'ID2', 'Value']!")
  
  SubjID_related = SparseGRMData %>% filter(ID1 != ID2) %>% select(ID1, ID2) %>%
    unlist(use.names = FALSE) %>% unique
  
  if(length(SubjID_related) == 0)
  {
    PairwiseIBD = as.data.frame(matrix(NA, 0, 5))
    
    colnames(PairwiseIBD) = c("ID1", "ID2", "pa", "pb", "pc")
    
    data.table::fwrite(PairwiseIBD, PairwiseIBDFile,
                       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }else
  {
    # read in the bim and fam data.
    bim = data.table::fread(bimFile)
    fam = data.table::fread(famFile)
    
    # check SNPs of frqFile and subjects of SparseGRMFile correspond to PlinkFile.
    if(any(GenoInfoMat$SNP != bim$V2))
      stop("Please check if SNPs in FrqFile match with PlinkFile")
    
    if(!all(SubjID_related %in% fam$V2))
      stop("Please check if related subjects in SparseGRMFile but not in PlinkFile")
    
    GenoInfoMat = GenoInfoMat %>% filter(MAF > minMafIBD)
    
    cat("Analyzing Number of Subjects:\t", length(SubjID_related), "\n")
    cat("Remaining Number of Markers:\t", nrow(GenoInfoMat), "\n")
    
    if(nrow(GenoInfoMat) < 1e4)
      warning("Number of Markers is a bit small, we recommend nSNPs > 10,000.\n")
    
    if(nrow(GenoInfoMat) > 1e5)
      warning("Number of Markers is a bit large, we recommend maxSampleNums < 2,500.\n")
    
    # metrics used in pairwise IBD calculation.
    altFreq = GenoInfoMat$MAF
    pro_var = 2 * (altFreq * (1 - altFreq))^2 # 2*pi^2*(1-pi)^2, where pi comes from Binom(2, pi).
    wi = sqrt(pro_var/(1 - pro_var)) # weights of each SNP
    
    # write the passed SNPIDs into IDsToInclude.
    IDsToIncludeFile = paste0(tempDir, "/IDsToInclude.txt")
    data.table::fwrite(data.table::data.table(GenoInfoMat$SNP), IDsToIncludeFile,
                       row.names = F, col.names = F, quote = F, sep = "\t")
    
    edges = t(SparseGRMData[, c("ID1", "ID2")])
    graph_GRM = make_graph(edges, directed = F)
    graph_list_all = graph_GRM %>% decompose()
    graph_length = lapply(graph_list_all, length)
    
    graph_list = graph_list_all[graph_length > 1]
    graph_length = lapply(graph_list, length) %>% unlist
    
    # initialize parameters
    PairwiseIBD = c();
    tSampleNums = 0;
    tSampleIDs = c();
    nParts = 1;
    
    # cycle for calculating pairwise IBD
    for(i in 1:length(graph_list))
    {
      tSampleNums = tSampleNums + graph_length[i];
      tSampleIDs = c(tSampleIDs, V(graph_list[[i]])$name);
      
      if(tSampleNums >= maxSampleNums | i == length(graph_list))
      {
        cat("\nProcessing the", nParts, "block(s) of Samples.\n")
        
        GenoList = GRAB.ReadGeno(bedFile,
                                 SampleIDs = tSampleIDs, 
                                 control = list(IDsToIncludeFile = IDsToIncludeFile,
                                                ImputeMethod = "mean"))
        
        tempGRM = SparseGRMData %>% 
          filter(ID1 %in% tSampleIDs & ID2 %in% tSampleIDs) %>%
          filter(ID1 != ID2) %>%
          mutate(idxID1 = match(ID1, rownames(GenoList$GenoMat))) %>%
          mutate(idxID2 = match(ID2, rownames(GenoList$GenoMat)))
        
        for(j in 1:nrow(tempGRM))
        {
          tempmetrics = GenoList$GenoMat[tempGRM$idxID1[j],] - GenoList$GenoMat[tempGRM$idxID2[j],]
          
          pc = 0.5 * weighted.mean(((abs(tempmetrics - 1) + abs(tempmetrics + 1)- 2)/pro_var), wi, na.rm = T)
          pc = ifelse(pc > (1-tempGRM$Value[j])^2, (1-tempGRM$Value[j])^2-1e-10, ifelse(pc < 1-2*tempGRM$Value[j], 1-2*tempGRM$Value[j], pc))
          
          pb = 2 - 2*pc - 2*tempGRM$Value[j]
          
          pa = 2*tempGRM$Value[j] + pc - 1
          
          if(pb < 0)
          {
            pa = pa + 0.5*pb; pb = 0; pc = 0;
          }
          
          PairwiseIBD = rbind(PairwiseIBD,
                              c(ID1 = tempGRM$ID1[j], ID2 = tempGRM$ID2[j], pa = pa, pb = pb, pc = pc))
        }
        
        cat("Completed analyzing the", nParts, "block(s) of Samples.\n")
        
        tSampleNums = 0; 
        tSampleIDs = c();
        nParts = nParts + 1;
      }
    }
    
    data.table::fwrite(data.table::data.table(PairwiseIBD), PairwiseIBDFile,
                       row.names = F, col.names = T, quote = F, sep = "\t")
    
    if(rm.tempFile)
      file.remove(IDsToIncludeFile)
  }
  
  message = paste("The PairwiseIBD has been stored in", PairwiseIBDFile)
  
  return(message)
}