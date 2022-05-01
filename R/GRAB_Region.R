
#' Conduct region-level genetic association testing
#' 
#' Test for association between phenotype of interest and regions including multiple genetic marker (mostly low-frequency or rare variants).
#' 
#' @param objNull the output object of function \code{\link{GRAB.NullModel}}. 
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the prefix is the same as GenoFile. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param OutputFile a character of output file to save the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point. If the program ends unexpectedly, the end point can help \code{GRAB} package understand where to restart the analysis. If \code{NULL} (default), \code{OutputFileIndex = paste0(OutputFile, ".index")}. 
#' @param RegionFile a character of region file to specify region-marker mapping with annotation information. Columns are separated by 'tab'. Column 1: region ID, Column 2: marker ID, Columns 3-n annotation (non-negative) similar as in STAAR. The header is required and the first two should be "REGION" and "MARKER".
#' @param RegionAnnoHeader a character vector of annotation in analysis. Optional, if not specified, all annotation columns are used in analysis. 
#' @param SparseGRMFile a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}.
#' @param MaxMAFVec a numeric vector of max MAF cutoff to include markers for region-level analysis. Default setting is \code{c(0.05, 0.01, 0.005)}.
#' @param control a list of parameters for controlling function \code{GRAB.Region}, more details can be seen in \code{Details} section.
#' @details 
#' \code{GRAB} package supports \code{SAIGE}, \code{POLMM}, and \code{SPACox} methods. 
#' Detailed information about the analysis methods is given in the \code{Details} section of \code{\link{GRAB.NullModel}}. 
#' Users do not need to specify them since functions \code{\link{GRAB.Marker}} and \code{GRAB.Region} will check the \code{class(objNull)}.
#' 
#' ## Region-based approaches are mostly for low-frequency and rare variants. 
#' 
#' ## The following details are about argument \code{control}
#' For PLINK files, the default \code{control$AlleleOrder = "alt-first"}; for BGEN files, the default \code{control$AlleleOrder = "ref-first"}.
#'   \itemize{
#'   \item \code{AlleleOrder}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   }
#'  The below is to customize the quality-control (QC) process.
#'   \itemize{
#'   \item \code{omp_num_threads}: (To be added later) a numeric value (default: value from data.table::getDTthreads()) to specify the number of threads in OpenMP for parallel computation.
#'   \item \code{ImputeMethod}: a character, "mean", "bestguess" (default), or "drop" (to be added later). Please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{MissingRateCutoff}: a numeric value *(default=0.15)*. Markers with missing rate > this value will be excluded from analysis.  
#'   \item \code{MinMACCutoff}: a numeric value *(default=5)*. Markers with MAC < this value will be treated as Ultra-Rare Variants (URV) and collapsed as one value.  
#'   \item \code{nRegionsEachChunk}: number of regions *(default=1)* in one chunk to output.
#'   }
#'   The below is for kernel-based approaches including SKAT and SKAT-O. For more details, please refer to \code{\link{SKAT}}.
#'   \itemize{
#'   \item \code{kernel}: a type of kernel *(default="linear.weighted")*.
#'   \item \code{weights_beta}: a numeric vector of parameters for the beta weights for the weighted kernels *(default=c(1, 25))*. 
#'   If you want to use your own weights, please use the \code{control$weights} parameter. It will be ignored if \code{control$weights} parameter is not \code{NULL}. 
#'   \item \code{weights}: a numeric vector of weights for the weighted kernels. If it is \code{NULL} (default), the beta weight with the \code{control$weights.beta} parameter is used.
#'   \item \code{r.corr}: the rho parameter for the compound symmetric correlation structure kernels. If you give a vector value, SKAT will conduct the optimal test. 
#'   It will be ignored if method="optimal" or method="optimal.adj" *(default=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1))*.  
#'   }
#'  The below is to customize the columns in the \code{OutputMarkerFile}. 
#'  Columns of \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, \code{MissingRate}, \code{Pvalue} are included for all methods.
#'  \itemize{
#'  \item \code{outputColumns}: For example, for POLMM method, users can set \code{control$outputColumns = c("beta", "seBeta", "AltFreqInGroup")}. 
#'     \itemize{
#'     \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; Optional: \code{zScore}, \code{AltFreqInGroup}, \code{nSamplesInGroup}, \code{AltCountsInGroup}
#'     \item \code{SPACox}: Optional: \code{zScore}
#'     }
#'  } 
#' @return Region-based analysis results are saved into two files: \code{OutputFile} and \code{OutputMarkerFile = paste0(OutputFile, ".markerInfo")}. 
#' 
#' The file of \code{OutputMarkerFile} is the same as the results of \code{\link{GRAB.Marker}}. The file of \code{OutputFile} includes columns as below.
#' \item{Region}{Region IDs from \code{RegionFile}}
#' \item{Anno.Type}{Annotation type from \code{RegionFile}}
#' \item{maxMAF}{the maximal cutoff of the MAF to select low-frequency/rare variants into analysis.}
#' \item{nSamples}{Number of samples in analysis.}
#' \item{nMarkers}{Number of markers whose MAF < \code{control$MaxMAFCutoff} and MAC > \code{control$MinMACCutoff}. Markers with annotation value <= 0 will be excluded from analysis.}
#' \item{nMarkersURV}{Number of Ultra-Rare Variants (URV) whose MAC < \code{control$MinMACCutoff}. Markers with annotation value <= 0 will be excluded from analysis.}
#' \item{pval.SKATO}{p-values based on SKAT-O method}
#' \item{pval.SKAT}{p-values based on SKAT method}
#' \item{pval.Burden}{p-values based on Burden test}
#' @examples 
#' objNullFile = system.file("results", "objNull.RData", package = "GRAB")
#' load(objNullFile)
#' class(objNull)    # "POLMM_NULL_Model", that indicates an object from POLMM method.
#' 
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/simuRegionOutput.txt")
#' GenoFile = system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
#' RegionFile = system.file("extdata", "simuRegion.txt", package = "GRAB")
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' 
#' ## make sure the output files does not exist at first
#' file.remove(OutputFile)
#' file.remove(paste0(OutputFile, ".markerInfo"))
#' file.remove(paste0(OutputFile, ".index"))
#' 
#' GRAB.Region(objNull, 
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile,
#'             RegionFile = RegionFile,
#'             SparseGRMFile = SparseGRMFile)
#'             
#' data.table::fread(OutputFile)
#' 
#' ## additional columns of "zScore", "nSamplesInGroup", "AltCountsInGroup", "AltFreqInGroup"
#' ## We do not recommend adding too many columns for all markers
#' 
#' file.remove(OutputFile)
#' file.remove(paste0(OutputFile, ".markerInfo"))
#' file.remove(paste0(OutputFile, ".index"))
#' GRAB.Region(objNull, 
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile,
#'             RegionFile = RegionFile,
#'             RegionAnnoHeader = c("ANNO1", "ANNO2"),
#'             SparseGRMFile = SparseGRMFile,
#'             control = list(outputColumns = c("beta", "seBeta", "zScore","nSamplesInGroup","AltCountsInGroup","AltFreqInGroup")))
#'             
#' data.table::fread(OutputFile)
#'     
#' @export
#' @import SKAT, data.table

GRAB.Region = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,
                       GroupFile,              # column 1: marker Set ID, column 2: SNP ID, columns 3-n: Annotations similar as in STAAR
                       # RegionAnnoHeader = NULL,
                       SparseGRMFile = NULL,
                       MaxMAFVec = c(0.05, 0.01, 0.005),
                       annoVec = c("stoploss:nonframeshift insertion:nonsynonymous SNV",
                                   "stoploss:nonframeshift insertion:nonsynonymous SNV:synonymous SNV:frameshift deletion:nonframeshift deletion:frameshift insertion:stopgain:splicing:startloss"),
                       chrom = "LOCO=F",
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);  # Check "Util.R"
  method = gsub("_NULL_Model", "", NullModelClass)
  
  if(is.null(OutputFileIndex)) 
    OutputFileIndex = paste0(OutputFile, ".index")
  
  outList = checkOutputFile(OutputFile, OutputFileIndex, "Region", 
                            nEachChunk = 1) # Check 'Util.R'
  
  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End
  
  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results have been saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    cat(message)
    return(message)
  }
  
  if(!Start){
    message = paste0("Parts of analysis have been conducted based on the index file:\n",
                     OutputFileIndex,"\n",
                     "The analysis will be restarted from chunk:\t",indexChunk+1,"\n");
    cat(message)
  }
  
  ## Check "control.R": if the setting of control is not specified, the default setting will be used
  control = checkControl.Region(control, NullModelClass)
  
  MaxMAF = max(MaxMAFVec)
  if(MaxMAF > 0.05) 
    stop("Maximal value of 'MaxMAFVec' should be <= 0.05.")
  control$max_maf_region = MaxMAF
  
  annoList = annoVec %>% 
    strsplit(split = ":") 
  
  allAnno = annoList %>% 
    unlist() %>% 
    unique()
  
  subjData = as.character(objNull$subjData);
  n = length(subjData)
  
  # note on 2022-04-26: check it later
  # Group = makeGroup(objNull$yVec)
  # ifOutGroup = any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)


  ## set up an object for genotype data
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData, control)  # Check 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  
  ## annotation in region
  # RegionList = getRegionList(RegionFile, RegionAnnoHeader, markerInfo)
  
  RegionList = getInfoGroupFile(GroupFile)
  nRegions = length(RegionList)
  
  # P1Mat = matrix(0, control$max_markers_region, n);
  # P2Mat = matrix(0, n, control$max_markers_region);
  
  # chrom1 = "FakeCHR";
  for(i in (indexChunk+1):nRegions){
    
    region = RegionList[[i]]
    
    regionID = region$regionID
    regionInfo = region$regionInfo
    
    regionInfo = markerInfo %>% 
      select(ID, genoIndex) %>% 
      merge(regionInfo, by = "ID") %>% 
      arrange(genoIndex) %>%
      filter(Annos %in% allAnno)
    
    nMarkers = nrow(regionInfo)
    if(nMarkers == 0)
      stop("nrow(regionInfo) == 0: no markers are found for region '", regionID, "'.")
    
    nAnno = length(annoList)
    annoMat = matrix(0, nrow = nMarkers, ncol = nAnno)
    colnames(annoMat) = annoVec
    
    for(iAnno in 1:nAnno)
      annoMat[,iAnno] = ifelse(regionInfo$Annos %in% annoList[[iAnno]], 1, 0)
    
    genoIndex = regionInfo$genoIndex
    weightVec = regionInfo$Weights
    
    weightExists = T  # update weight parts later: 2022-05-01
    if(all(is.na(weightVec))){
      weightExists = F
      weightVec = rep(1, nMarkers)
    }else{
      if(any(is.na(weightVec) | weightVec <= 0))
        stop("The provided marker weights cannot be non-positive (<= 0) or NA.")
    }
    
    print(paste0("Analyzing Region of ", regionID, " (",i,"/",nRegions,")."))
    print(paste(regionInfo$ID, collapse = ", "))
    
    with(control, 
         setRegion_GlobalVarsInCPP(impute_method, 
                                   missing_cutoff, 
                                   max_maf_region, 
                                   min_mac_region, 
                                   max_markers_region, 
                                   omp_num_threads))
    textToParse = paste0("obj.setRegion = setRegion.", method, "(objNull, control, chrom, SparseGRMFile)")
    eval(parse(text = textToParse))
    
    subjLabel = rep(0, n)  # 2022-04-27: give labels to each subject (e.g. 0 for control and 1 for case), to be extended later. Start from 0.
    nLabel = max(subjLabel) + 1
    
    obj.mainRegionInCPP = mainRegionInCPPcheck(method, genoType, genoIndex, weightVec, OutputFile, 
                                               subjLabel, nLabel, 
                                               annoMat, annoVec)
    
    textToParse = paste0("obj.mainRegion = mainRegion.", method, ".Check(genoType, genoIndex, OutputFile, control, n, obj.setRegion, obj.mainRegionInCPP)")
    eval(parse(text = textToParse))
    
    Other.Markers = obj.mainRegion$Other.Markers %>% mutate(Region = regionID, .before = ID)
    VarMat = obj.mainRegion$VarMat
    RV.Markers0 = obj.mainRegion$RV.Markers %>% mutate(Region = regionID, .before = ID)
    
    if(nrow(VarMat) != nrow(RV.Markers0))
      stop("nrow(VarMat) != nrow(RV.Markers0)!")
    
    RV.Markers = RV.Markers0 %>% 
      mutate(betaWeights = dbeta(MAF, control$weights.beta[1], control$weights.beta[2]),
             adjVarSVec = StatVec^2 / qchisq(pval1Vec, df = 1, lower.tail = F),
             r0 = adjVarSVec / diag(VarMat),
             wr0 = sqrt(r0) * betaWeights,
             wStatVec = StatVec * betaWeights)
    # check given weights version later: 2022-05-01
    
    wr0 = RV.Markers$wr0
    
    wadjVarSMat = t(VarMat * wr0) * wr0
    
    # RV.Markers %>% head()
    
    RV.MarkersWithAnno = regionInfo %>% 
      select(-genoIndex) %>% 
      merge(RV.Markers %>% select(ID, MAF, posRow), by = "ID")
    
    Other.MarkersWithAnno = regionInfo %>% 
      select(ID, Annos) %>% 
      merge(Other.Markers %>% filter(IndicatorVec == 2) %>% select(ID), by = "ID")
    
    RV.MarkersURV = RV.Markers %>% filter(Info == "Ultra-Rare Variants") %>% select(ID, posRow)
    
    pval.Region = data.frame()
    for(anno in annoVec)
    {
      posURV = RV.MarkersURV %>% filter(ID == anno) %>% select(posRow) %>% unlist()
      nMarkersURV = Other.MarkersWithAnno %>% filter(Annos %in% annoTemp) %>% nrow()
      if(length(posURV) != 1)
        stop("length(posURV) != 1")
      
      annoTemp = unlist(strsplit(anno, split = ":"))
      for(MaxMAF in MaxMAFVec)
      {
        posRV = RV.MarkersWithAnno %>% filter(MAF < MaxMAF & Annos %in% annoTemp) %>% select(posRow) %>% unlist()
        pos = c(posRV, posURV)
        n1 = length(pos)
        out_SKAT_List = with(RV.Markers, try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec[pos], 
                                                                        Phi = wadjVarSMat[pos, pos],  
                                                                        r.corr = control$r.corr, 
                                                                        method = "optimal.adj", 
                                                                        Score.Resampling = NULL),
                                             silent = TRUE))
        
        if(class(out_SKAT_List) == "try-error"){
          Pvalue = c(NA, NA, NA)
          error.code = 2
        }else if(!any(c(0,1) %in% out_SKAT_List$param$rho)){
          Pvalue = c(NA, NA, NA)
          error.code = 3
        }else{
          pos00 = which(out_SKAT_List$param$rho == 0)
          pos01 = which(out_SKAT_List$param$rho == 1)
          Pvalue = c(out_SKAT_List$p.value,                  # SKAT-O
                     out_SKAT_List$param$p.val.each[pos00],   # SKAT
                     out_SKAT_List$param$p.val.each[pos01])   # Burden Test
          error.code = 0
        }
        
        pval.Region = rbind.data.frame(pval.Region,
                                       data.frame(Region = regionID,
                                                  nMarkers = length(posRV),
                                                  nMarkersURV = nMarkersURV,
                                                  Anno.Type = anno,
                                                  MaxMAF.Cutoff = MaxMAF,
                                                  pval.SKATO = Pvalue[1], 
                                                  pval.SKAT = Pvalue[2],
                                                  pval.Burden = Pvalue[3]))
      }
    }
    
    writeOutputFile(Output = list(pval.Region, 
                                  RV.Markers0, 
                                  Other.Markers), 
                    OutputFile = list(OutputFile, 
                                      paste0(OutputFile, ".markerInfo"), 
                                      paste0(OutputFile, ".otherMarkerInfo")), 
                    OutputFileIndex = OutputFileIndex,
                    AnalysisType = "Region",
                    nEachChunk = 1,
                    indexChunk = i,
                    Start = (i==1),
                    End = (i==nRegions))
  }
      
  cat("Analysis done! The results have been saved to the below files:\n", 
      OutputFile, "\n",
      paste0(OutputFile, ".markerInfo\n"),
      paste0(OutputFile, ".otherMarkerInfo\n"))
  
  message = paste0("Analysis done! The main results have been saved to '", OutputFile,"'")
  
  return(message)
}


# setRegion = function(NullModelClass, objNull, control, chrom, SparseGRMFile, Group, ifOutGroup)
# {
#   # The following function is in Main.cpp
#   nGroup = length(unique(Group))
#   setRegion_GlobalVarsInCPP(control$impute_method,
#                             control$missing_cutoff,
#                             control$max_maf_region,
#                             control$min_mac_region,
#                             control$max_markers_region,
#                             control$omp_num_threads,
#                             Group, ifOutGroup, nGroup)
#   
#   # Check POLMM.R
#   if(NullModelClass == "POLMM_NULL_Model")
#     obj.setRegion = setRegion.POLMM(objNull, control, chrom, SparseGRMFile)  
#   
#   # To be continued
#   if(NullModelClass == "SAIGE_NULL_Model")
#     obj.setRegion = setRegion.SAIGE(objNull, control, chrom, SparseGRMFile)
#   
#   # Check SPACox.R
#   if(NullModelClass == "SPACox_NULL_Model")
#     obj.setRegion = setRegion.SPACox(objNull, control)
#   
#   return(obj.setRegion)
# }

# mainRegion = function(NullModelClass, genoType, genoIndex, OutputFile, n, P1Mat, P2Mat, outputColumns)
# {
#   method = gsub("_NULL_Model$", "", NullModelClass)
#   
#   obj.mainRegion = mainRegionInCPP(method, genoType, genoIndex, OutputFile, n, P1Mat, P2Mat)
#   
#   ## required columns for all methods
#   info.Region = with(obj.mainRegion, data.frame(Marker = markerVec,
#                                                 Info = infoVec,
#                                                 AltFreq = altFreqVec,
#                                                 MAC = MACVec,
#                                                 MAF = MAFVec,
#                                                 MissingRate = missingRateVec, 
#                                                 IsUltraRareVariants = indicatorVec - 1,
#                                                 stringsAsFactors = F))
#   
#   if(NullModelClass == "POLMM_NULL_Model")
#   {
#     optionalColumns = c("beta", "seBeta", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
#     additionalColumns = intersect(optionalColumns, outputColumns)
#     
#     if(length(additionalColumns) > 0)
#       info.Region = cbind.data.frame(info.Region, 
#                                      as.data.frame(obj.mainRegion[additionalColumns]))
#   }
#     
#   if(NullModelClass == "SPACox_NULL_Model")
#   {
#     # obj.mainRegion = mainRegionInCPP("SPACox", genoType, genoIndex, OutputFile, n, P1Mat, P2Mat)
#   }
#   
#   ### remove rows whose markers do not pass QC
#   
#   info.Region = subset(info.Region, IsUltraRareVariants != -1)
#   pos = which(info.Region$IsUltraRareVariants == 0)
#   
#   info.Region$Pvalue = NA
#   info.Region$Pvalue[pos] = obj.mainRegion$pval1Vec
#   
#   return(list(StatVec = obj.mainRegion$StatVec,
#               pval1Vec = obj.mainRegion$pval1Vec,
#               VarMat = obj.mainRegion$VarMat,
#               info.Region = info.Region))
# }
# 
# mainRegionURV = function(NullModelClass,
#                          genoType,
#                          genoIndex,
#                          n)
# {
#   if(NullModelClass == "POLMM_NULL_Model")
#     obj.mainRegionURV = mainRegionURVInCPP("POLMM", genoType, genoIndex, n)
#   
#   if(NullModelClass == "SPACox_NULL_Model")
#     obj.mainRegionURV = mainRegionURVInCPP("SPACox", genoType, genoIndex, n)
#   
#   return(obj.mainRegionURV)
# }

# Rcpp::List mainRegionInCPP(std::string t_method,       // "POLMM", "SAIGE"
#                            std::string t_genoType,     // "PLINK", "BGEN"
#                            std::vector<uint32_t> t_genoIndex,
#                            std::string t_outputFile,
#                            unsigned int t_n)           // sample size  

# mainRegion = function(NullModelClass, genoType, genoIndex, regionMat)
# {
#   # the following function is in POLMM.R
#   if(NullModelClass == "POLMM_NULL_Model")
#     obj.mainRegion = mainRegion.POLMM(genoType, genoIndex, regionMat)
#   
#   if(NullModelClass == "SAIGE_NULL_Model")
#     obj.mainRegion = mainRegion.SAIGE(genoType, genoIndex, regionMat)
#   
#   if(NullModelClass == "SPACox_NULL_Model")
#     obj.mainRegion = mainRegion.SPACox(genoType, genoIndex, regionMat)
#   
#   return(obj.mainRegion)
# }

# setRegion = function()
# {
#   if(objNull$controlList$LOCO){
#     if(!chrom %in% names(objNull$LOCOList))
#       stop("'chrom' should be in names(objNull$LOCOList).")
#     obj.CHR = objNull$LOCOList[[chrom]]
#     
#     if(!chrom %in% names(SparseGRM))
#       stop("'chrom' should be in names(SparseGRM).")
#     SparseGRM.CHR = SparseGRM[[chrom]]
#   }
#   
#   SPmatR.CHR = makeSPmatR(SparseGRM.CHR, objNull$subjIDs)
#   
#   setPOLMMobjInR(obj.CHR$muMat,
#                  obj.CHR$iRMat,
#                  objNull$Cova,
#                  objNull$yVec,          # 1 to J
#                  SPmatR.CHR,
#                  objNull$tau,
#                  POLMM.control$printPCGInfo,
#                  POLMM.control$tolPCG,
#                  POLMM.control$maxiterPCG)
#   
#   memory_chunk = POLMM.control$memory_chunk
#   
#   # to be continued
#   n = length(SubjID.step1)
#   J = max(objNull$yVec)
#   p = ncol(objNull$Cova)
#   NonZero_cutoff = floor(log(1e7, J))  # for efficient resampling (ER)
#   maxMarkers = getMaxMarkers(memory_chunk, n, J, p);
#   
#   print(paste0("The current POLMM.control$memory_chunk is ", memory_chunk,"(GB)."))
#   print(paste0("Based on the sample size, we divide region with more than ", maxMarkers, " markers into multiple chunks to save memory usage."))
#   print("If the memory usage still exceed the memory you request, please set a smaller POLMM.control$memory_chunk.")
#   
#   StdStat_cutoff = POLMM.control$SPA_cutoff;
# }

# calculate region-based p-values for any given annotations
# mainRegion = function(region, NullModelClass, objNull, control)
# {
#   markers = region$markers
#   annoMat = region$annoMat
#   
#   OutList = mainRegioninCPP(markers, NullModelClass)
#   
#   ## extract information from mainRegioninCPP(.)
#   StatVec = OutList$StatVec
#   VarSVec = diag(OutList$VarSMat)
#   pvalNormVec = OutList$pvalNormVec;
#   adjPVec = OutList$pvalVec;
#   adjVarSVec = StatVec^2 / qchisq(adjPVec, df = 1, lower.tail = F)
#   weightVec = OutList$weightVec
#   
#   # index of SNPs passing criterion (MAF, missing rate, et al.)
#   posVec = OutList$posVec + 1   # "+1" because C++ starts from 0 and R starts from 1
#   annoMat = annoMat[posVec,,drop=F]
#   q = ncol(annoMat)             # number of annotations: column 1 is always 1s
#   
#   r0 = adjVarSVec / VarSVec 
#   r0 = pmax(r0, 1)
#   
#   pvalVec.BT = c()
#   pvalVec.SKAT = c()
#   pvalVec.SKATO = c()
#   annoVec = c()
#   errVec = c()
#   
#   # cycle for multiple annotations
#   for(i in 1:q){
#     annoName = colnames(annoMat)[i]
#     annoWeights = weightVec * annoMat[,i]
#     wr0 = sqrt(r0) * annoWeights
#     wStatVec = StatVec * annoWeights
#     wadjVarSMat = t(OutList$VarSMat * wr0) * wr0
#     wadjVarSMat = wadjVarSMat * max(OutList$rBT, 1)
#     
#     out_SKAT_List = try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec, 
#                                                    Phi = wadjVarSMat,
#                                                    r.corr = control$r_corr, 
#                                                    method = "optimal.adj", 
#                                                    Score.Resampling = NULL),
#                         silent = TRUE)
#     
#     # betaVec = StatVec / adjVarSVec; 
#     if(class(out_SKAT_List) == "try-error"){
#       Pvalue = c(NA, NA, NA)
#       errCode = 2
#     }else if(!any(c(0,1) %in% out_SKAT_List$param$rho)){
#       Pvalue = c(NA, NA, NA)
#       errCode = 3
#     }else{
#       pos0 = which(out_SKAT_List$param$rho == 0)
#       pos1 = which(out_SKAT_List$param$rho == 1)
#       Pvalue = c(out_SKAT_List$p.value,                  # SKAT-O
#                  out_SKAT_List$param$p.val.each[pos0],   # SKAT
#                  out_SKAT_List$param$p.val.each[pos1])   # Burden Test
#       errCode = 0
#     }
#     
#     pvalVec.SKATO = c(pvalVec.SKATO, Pvalue[1])
#     pvalVec.SKAT = c(pvalVec.SKAT, Pvalue[2])
#     pvalVec.BT = c(pvalVec.BT, Pvalue[3])
#     annoVec = c(annoVec, annoName)
#     errVec = c(errVec, errCode)
#   }
#   
#   ###################
#   
#   OUT.Region = rbind(OUT.Region,
#                      c(Region, 
#                        length(OutList$markerVec), 
#                        paste(annoVec, collapse = ","),
#                        paste(pvalVec.SKATO, collapse = ","),
#                        paste(pvalVec.SKAT, collapse = ","),
#                        paste(pvalVec.BT, collapse = ","),
#                        paste(errVec, collapse = ","),
#                        paste(OutList$markerVec, collapse = ","),
#                        paste(OutList$freqVec, collapse = ","),
#                        paste(OutList$flipVec, collapse = ","),
#                        paste(StatVec, collapse = ","),
#                        paste(adjVarSVec, collapse = ","),
#                        paste(adjPVec, collapse = ",")))
#   
#   tempFiles = list.files(path = dirname(OutputFile),
#                          pattern = paste0("^",basename(OutputFile),".*\\.bin$"),
#                          full.names = T)
#   file.remove(tempFiles)
#   
#   # out_Multi_Set = rnorm(1)
#   colnames(OUT.Region) = c("regionName", "nMarkers", "Annotation", "P.SKAT-O", "P.SKAT", "P.Burden",
#                            "error.code", "markerInfo", "markerMAF","markerAlleleFlip",
#                            "markerStat","markerVarS","markerPvalue")
#   
#   return(resRegion)
# }


# extract region-marker mapping from regionFile
# getRegionList = function(RegionFile,
#                          RegionAnnoHeader,
#                          markerInfo)
# {
#   cat("Start extracting marker-level information from 'RegionFile' of", RegionFile, "....\n")
#   
#   if(!file.exists(RegionFile))
#     stop("Cannot find 'RegionFile' in ", RegionFile)
#   
#   RegionData = data.table::fread(RegionFile, header = T, stringsAsFactors = F, sep = "\t");
#   RegionData = as.data.frame(RegionData)
#   colnames(RegionData)[1:2] = toupper(colnames(RegionData)[1:2])
#   
#   if(any(colnames(RegionData)[1:2] != c("REGION", "MARKER")))
#     stop("The first two elements in the header of 'RegionFile' should be c('REGION', 'MARKER').")
#   
#   # updated on 2021-08-05
#   colnames(markerInfo)[3] = "MARKER"
#   RegionData = merge(RegionData, markerInfo, by = "MARKER", all.x = T, sort = F)
#   posNA = which(is.na(RegionData$genoIndex))
#   
#   if(length(posNA) != 0){
#     print(head(RegionData[posNA,1:2]))
#     stop("Total ",length(posNA)," markers in 'RegionFile' are not in 'GenoFile'. 
#          Please remove these markers before region-level analysis.")
#   }
#   
#   HeaderInRegionData = colnames(RegionData)
#   if(!is.null(RegionAnnoHeader)){
#     if(any(!RegionAnnoHeader %in% HeaderInRegionData))
#       stop("At least one element in 'RegionAnnoHeader' is not in the header of RegionFile.")
#     posAnno = which(HeaderInRegionData %in% RegionAnnoHeader)
#   }else{
#     print("Since no 'RegionAnnoHeader' is given, region-based testing will not incorporate any annotation information.")
#     posAnno = NULL
#   }
#   
#   RegionList = list()
#   uRegion = unique(RegionData$REGION)
#   for(r in uRegion){
#     
#     # print(paste0("Analyzing region ",r,"...."))
#     
#     posSNP = which(RegionData$REGION == r)
#     SNP = RegionData$MARKER[posSNP]
#     
#     if(any(duplicated(SNP)))
#       stop("Please check RegionFile: in region ", r,": duplicated SNPs exist.")
#     
#     # posMarker = match(SNP, markerInfo$ID, 0)
#     # if(any(posMarker == 0))
#     #   stop(paste0("At least one marker in region ", r," are not in 'GenoFile' and 'GenoFileIndex'."))
#     
#     regionMat = cbind(BASE=1, RegionData[posSNP, posAnno, drop=F])
#     rownames(regionMat) = SNP
#     
#     # genoIndex = markerInfo$genoIndex[posMarker]
#     # chrom = markerInfo$CHROM[posMarker]
#     genoIndex = RegionData$genoIndex[posSNP]
#     chrom = RegionData$CHROM[posSNP]
#     uchrom = unique(chrom)
#     
#     if(length(uchrom) != 1)
#       stop("In region ",r,", markers are from multiple chromosomes.")
#     
#     RegionList[[r]] = list(SNP = SNP,
#                            regionMat = regionMat,
#                            genoIndex = genoIndex,
#                            chrom = uchrom)
#   }
#   
#   return(RegionList)
# }

getInfoGroupLine = function(markerGroupLine, nLine)
{
  if(length(markerGroupLine) == 0)
    stop("The line ", nLine," in `groupFile` is empty.")
  
  info = unlist(strsplit(markerGroupLine, split = "\t"))
  if(length(info) < 3)
    stop("The line ", nLine, " in `groupFile` includes < 3 elements, please note that each line should be seperated by 'tab'.")
  
  geneID = info[1];
  type = info[2];
  values = info[c(-1,-2)]
  
  if(!type %in% c("var", "anno", "weight"))
    stop("The second column of the groupFile (tab-seperated) should be one of 'var', 'anno', and 'weight'.\n
         Please double check line ", nLine, ".")
  
  if(type == "weight")
    values = as.numeric(values)
  
  n = length(values)
  infoList = list(geneID = geneID,
                  type = type,
                  values = values,
                  n = n)
  return(infoList)
}

getInfoGroupFile = function(GroupFile)
{
  cat("Start extracting marker-level information from 'GroupFile':\n", 
      GroupFile, "\n")	
  if(!file.exists(GroupFile))
    stop("cannot find the below file:\n", GroupFile)
  
  gf = file(GroupFile, "r")
  nRegion = 0
  regionList = list()
  nLine = 1
  
  previousType = "first";
  previousGene = "first";
  Weights = NA
  nRegion = 1;
  
  while(TRUE)
  {
    markerGroupLine = readLines(gf, n = 1)
    
    if(length(markerGroupLine) == 0)
    {
      if(nRegion == 0)
        stop("Cannot find any region information in 'GroupFile'.")
      regionList[[nRegion]] = list(regionID = previousGene,
                                   regionInfo = data.frame(ID = Markers,
                                                           Annos = Annos,
                                                           Weights = Weights))
      close(gf);
      break;
    }
    
    infoList = getInfoGroupLine(markerGroupLine, nLine)
    nLine = nLine + 1
    
    geneID = infoList$geneID
    type = infoList$type
    values = infoList$values
    n = infoList$n
    
    if(type == "var")
    {
      if(previousType == "var")
        stop("Cannot find 'anno' line for region ", previousGene, ".")
      if(previousType != "first"){
        regionList[[nRegion]] = list(regionID = previousGene,
                                     regionInfo = data.frame(ID = Markers,
                                                             Annos = Annos,
                                                             Weights = Weights))
        nRegion = nRegion + 1;
      }
      
      # cat("Region ", geneID, " includes ", n, "variants.\n")
      Markers = values
      n1 = n
      Weights = NA
    }
    
    if(type == "anno")
    {
      if(n != n1)
        stop("The length of annotations for markers is not equal to the length of marker IDs")
      if(previousType != "var")
        stop("In the 'GroupFile', the 'anno' line should follow the 'var' line.")
      Annos = values;
    }
    
    if(type == "weight")
    {
      if(n != n1)
        stop("The length of weights for markers is not equal to the length of marker IDs")
      if(previousType != "anno")
        stop("In the 'GroupFile', the 'weight' line should follow the 'anno' line.")
      Weights = values
    }
    
    previousType = type
    previousGene = geneID
  }
  
  cat("Total ", nRegion, " groups are in 'GroupFile'.\n\n")
  return(regionList)
}
