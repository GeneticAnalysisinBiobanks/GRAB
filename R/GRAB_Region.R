
#' Conduct region-level genetic association testing
#' 
#' Test for association between phenotype of interest and regions including multiple genetic marker (mostly low-frequency or rare variants).
#' 
#' @param objNull the output object of function \code{\link{GRAB.NullModel}}. 
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the prefix is the same as GenoFile. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param OutputFile a character of output file to save the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point. If the program ends unexpectedly, the end point can help \code{GRAB} package understand where to restart the analysis. If \code{NULL} (default), \code{OutputFileIndex = paste0(OutputFile, ".index")}. 
#' @param RegionFile a character of region file to specify region-marker mapping with annotation information. Column 1: region ID, Column 2: marker ID, Columns 3-n annotation (non-negative) similar as in STAAR. The header is required and the first two should be "REGION" and "MARKER".
#' @param RegionAnnoHeader a character vector of annotation in analysis. Optional, if not specified, all annotation columns are used in analysis. 
#' @param control a list of parameters for controlling function \code{GRAB.Region}, more details can be seen in \code{Details} section.
#' @param SparseGRMFile a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}.
#' @details 
#' \code{GRAB} package supports \code{SAIGE}, \code{POLMM}, and \code{SPACox} methods. 
#' Detailed information about the analysis methods is given in the \code{Details} section of \code{\link{GRAB.NullModel}}. 
#' Users do not need to specify them since functions \code{\link{GRAB.Marker}} and \code{GRAB.Region} will check the \code{class(objNull)}.
#' 
#' Region-based approaches are mostly for low-frequency and rare variants. 
#' 
#' ## The following details are about argument \code{control}
#' For PLINK files, the default \code{control$AlleleOrder = "alt-first"}; for BGEN files, the default \code{control$AlleleOrder = "alt-first"}.
#'   \itemize{
#'   \item \code{AlleleOrder}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   }
#'  The below is to customize the quality-control (QC) process.
#'   \itemize{
#'   \item \code{omp_num_threads}: (To be added later) a numeric value (default: value from data.table::getDTthreads()) to specify the number of threads in OpenMP for parallel computation.
#'   \item \code{ImputeMethod}: a character, "mean", "bestguess", or "drop" (to be added later). Please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{MissingRateCutoff}: a numeric value *(default=0.15)*. Markers with missing rate > this value will be excluded from analysis.  
#'   \item \code{MaxMAFCutoff}: a numeric value *(default=0.01)*. Markers with MAF > this value will be excluded from analysis.  
#'   \item \code{MinMACCutoff}: a numeric value *(default=10)*. Markers with MAC < this value will be treated as Ultra-Rare Variants (URV) and collapsed as one value.  
#'   \item \code{nMarkersEachChunk}: number of markers *(default=300)* in one chunk to output.
#'   }
#'   The below is for kernel-based approaches including SKAT and SKAT-O. For more details, please refer to \code{\link{SKAT}}.
#'   \itemize{
#'   \item \code{kernel}: a type of kernel (default= "linear.weighted"). *(default="linear.weighted")*
#'   \item \code{weights_beta}: a numeric vector of parameters for the beta weights for the weighted kernels. If you want to use your own weights, please use the \code{control$weights} parameter. It will be ignored if \code{control$weights} parameter is not \code{NULL}. *(default=c(1,25))*
#'   \item \code{weights}: a numeric vector of weights for the weighted kernels. If it is \code{NULL}, the beta weight with the \code{control$weights.beta} parameter is used. *(default=NULL)*
#'   \item \code{r.corr}: the rho parameter for the compound symmetric correlation structure kernels. If you give a vector value, SKAT will conduct the optimal test. It will be ignored if method=“optimal” or method=“optimal.adj” *(default=c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1))*.  
#'   }
#'  The below is to customize the columns in the \code{OutputFile}. 
#'  Columns of \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, \code{MissingRate}, \code{Pvalue} are included for methods.
#'  \itemize{
#'  \item \code{outputColumns}: For example, for POLMM method, users can set \code{control$outputColumns = c("beta", "seBeta", "AltFreqInGroup")}. 
#'     \itemize{
#'     \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; Optional: \code{zScore}, \code{altFreqInGroup}
#'     \item \code{SPACox}: Optional: \code{PvalueNorm}, \code{zScore}
#'     }
#'  } 
#' @return Region-based analysis results are saved in two files: \code{OutputFile} and \code{OutputMarkerFile = paste0(OutputFile, ".markerInfo")}. 
#' 
#' The file of \code{OutputMarkerFile} is the same as the results of \code{\link{GRAB.Marker}}. The file of \code{OutputFile} includes columns as below.
#' \item{Region}{Region IDs from \code{RegionFile}}
#' \item{nMarkers}{Number of markers whose MAF < \code{control$MaxMAFCutoff} and MAC > \code{control$MinMACCutoff}.}
#' \item{nMarkersURV}{Number of Ultra-Rare Variants (URV) whose MAC < \code{control$MinMACCutoff}.}
#' \item{Anno.Type}{Annotation type from \code{RegionFile}}
#' \item{pval.SKATO}{p-values based on SKAT-O method}
#' \item{pval.SKAT}{p-values based on SKAT method}
#' \item{pval.Burden}{p-values based on Burden test}
#' @examples 
#' # We put examples to the specific help pages for different methods. 
#' # If you want to use "POLMM" method, please check ?GRAB.POLMM for more details.
#' @export
#' @import SKAT, data.table

GRAB.Region = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,
                       RegionFile,              # column 1: marker Set ID, column 2: SNP ID, columns 3-n: Annotations similar as in STAAR
                       RegionAnnoHeader = NULL,
                       SparseGRMFile = NULL,
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);  # this function is in "Util.R"
  
  if(is.null(OutputFileIndex)) 
    OutputFileIndex = paste0(OutputFile, ".index")
  
  outIndex = checkOutputFile(OutputFile, OutputFileIndex, "Region", 1) # this function is in 'Util.R'
  
  ## check the setting of control, if not specified, the default setting will be used
  control = checkControl.Region(control, NullModelClass)
  
  subjData = as.character(objNull$subjData);
  n = length(subjData)
  
  Group = makeGroup(objNull$yVec)
  if(any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)){
    ifOutGroup = TRUE
  }else{
    ifOutGroup = FALSE
  }
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData, control)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  
  ## annotation
  RegionList = getRegionList(RegionFile, RegionAnnoHeader, markerInfo)
  nRegions = length(RegionList)
  
  if(outIndex == nRegions + 1)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results can be seen in '", OutputFile, "'.",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }
    
  P1Mat = matrix(0, control$max_markers_region, n);
  P2Mat = matrix(0, n, control$max_markers_region);
  
  chrom1 = "FakeCHR";
  for(i in outIndex:nRegions){
    
    region = RegionList[[i]]
    regionName = names(RegionList)[i]
    
    SNP = region$SNP
    regionMat = region$regionMat  # annotation values
    genoIndex = region$genoIndex
    chrom = region$chrom
    
    print(paste0("Analyzing Region of ", regionName, " (",i,"/",nRegions,")."))
    print(paste(SNP, collapse = ", "))
    
    if(chrom1 != chrom){
      obj.setRegion = setRegion(NullModelClass, objNull, control, chrom, SparseGRMFile, Group, ifOutGroup)
      chrom1 = chrom
    }
    
    # main function to calculate summary statistics for region-based analysis 
    obj.mainRegion = mainRegion(NullModelClass, genoType, genoIndex, OutputFile, n, P1Mat, P2Mat)
    
    if(length(obj.mainRegion) == 0){
      writeOutputFile(list(), i, 
                      list(), 
                      OutputFileIndex, "Region", 1)
      next;
    }
    
    ###
    # MAF = pmin(obj.mainRegion$altFreqVec, 1-obj.mainRegion$altFreqVec)
    MAF = pmin(obj.mainRegion$MAFVec)
    weights = dbeta(MAF, control$weights.beta[1], control$weights.beta[2])
    
    StatVec = obj.mainRegion$StatVec
    
    VarSVec = diag(obj.mainRegion$VarMat)
    adjPVec = obj.mainRegion$adjPVec;
    adjVarSVec = StatVec^2 / qchisq(adjPVec, df = 1, lower.tail = F)
    
    r0 = adjVarSVec / VarSVec 
    r0 = pmax(r0, 1)
    
    # info.Region = data.frame(Region = regionName,
    #                          nMarker = length(obj.mainRegion$markerVec),
    #                          Markers = paste0(obj.mainRegion$markerVec, collapse = ","),
    #                          Info = paste0(obj.mainRegion$infoVec, collapse = ","),
    #                          AltFreq = paste0(obj.mainRegion$altFreqVec, collapse = ","),
    #                          MissingRate = paste0(obj.mainRegion$missingRateVec, collapse = ","),
    #                          # Stat = paste0(obj.mainRegion$StatVec, collapse = ","),
    #                          Beta = paste0(obj.mainRegion$BetaVec, collapse = ","),
    #                          seBeta = paste0(obj.mainRegion$seBetaVec, collapse = ","),
    #                          pval0 = paste0(obj.mainRegion$pval0Vec, collapse = ","),
    #                          pval1 = paste0(obj.mainRegion$pval1Vec, collapse = ","))
    
    nMarker = length(obj.mainRegion$markerVec)
    nMarkerURV = length(obj.mainRegion$markerURVVec)
    
    if(nMarker <= control$min_nMarker){
      writeOutputFile(list(), i, 
                      list(), 
                      OutputFileIndex, "Region", 1)
      next;
    }
    
    info.Marker.Region = data.frame(Region = regionName,
                                    Marker = obj.mainRegion$markerVec,
                                    IsUltraRareVariants = 0,
                                    Info = obj.mainRegion$infoVec,
                                    AltFreq = obj.mainRegion$altFreqVec,
                                    MAC = obj.mainRegion$MACVec,
                                    MissingRate = obj.mainRegion$missingRateVec,
                                    Beta = obj.mainRegion$BetaVec,
                                    seBeta = obj.mainRegion$seBetaVec,
                                    pval0 = obj.mainRegion$pval0Vec,
                                    pval1 = obj.mainRegion$pval1Vec)
    
    if(length(obj.mainRegion$markerURVVec) == 0){
      info.MarkerURV.Region = NULL;
    }else{
      info.MarkerURV.Region = data.frame(Region = regionName,
                                         Marker = obj.mainRegion$markerURVVec,
                                         IsUltraRareVariants = 1,
                                         Info = obj.mainRegion$infoURVVec,
                                         AltFreq = obj.mainRegion$altFreqURVVec,
                                         MAC = obj.mainRegion$MACURVVec,
                                         MissingRate = obj.mainRegion$missingRateURVVec,
                                         Beta = NA,
                                         seBeta = NA,
                                         pval0 = NA,
                                         pval1 = NA)
    }

    
    info.Region = rbind.data.frame(info.Marker.Region, info.MarkerURV.Region)
    
    posMarker = match(obj.mainRegion$markerVec, SNP)
    posMarkerURV = match(obj.mainRegion$markerURVVec, SNP)
    
    
    # print("length(posMarker)/length(posMarkerURV)/length(weights)/length(r0):")
    # print(c(length(posMarker), length(posMarkerURV), length(weights), length(r0)))
    # 
    # print("dim(obj.mainRegion$VarMat):")
    # print(dim(obj.mainRegion$VarMat))
    
    # print(obj.mainRegion$markerVec)
    # print(SNP)
    # print(posMarker)
    
    pval.Region = data.frame()
    for(j in 1:ncol(regionMat))
    {
      AnnoName = colnames(regionMat)[j]
      
      regionDataTemp1 = regionMat[posMarker, j]
      regionDataTemp2 = regionMat[posMarkerURV, j]
      regionData = c(regionDataTemp1, mean(regionDataTemp2))
      
      # print("length(regionDataTemp1)/length(regionDataTemp2)/length(regionData):")
      # print(c(length(regionDataTemp1), length(regionDataTemp2), length(regionData)))

      AnnoWeights = weights * regionData 
      
      # print(AnnoWeights)
      
      wr0 = sqrt(r0) * AnnoWeights
      
      # print("length(wr0)/length(AnnoWeights):")
      # print(c(length(wr0), length(AnnoWeights)))
      
      wStatVec = StatVec * AnnoWeights
      wadjVarSMat = t(obj.mainRegion$VarMat * wr0) * wr0
      
      # print("length(wStatVec):")
      # print(length(wStatVec))
      # 
      # print("dim(wadjVarSMat):")
      # print(dim(wadjVarSMat))
      
      # print(r0)
      # print(StatVec)
      # print(wStatVec)
      # print(wadjVarSMat)
      
      out_SKAT_List = try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec, 
                                                     Phi = wadjVarSMat,
                                                     r.corr = control$r.corr, 
                                                     method = "optimal.adj", 
                                                     Score.Resampling = NULL),
                          silent = TRUE)
      
      if(class(out_SKAT_List) == "try-error"){
        Pvalue = c(NA, NA, NA)
        error.code = 2
      }else if(!any(c(0,1) %in% out_SKAT_List$param$rho)){
        Pvalue = c(NA, NA, NA)
        error.code = 3
      }else{
        pos0 = which(out_SKAT_List$param$rho == 0)
        pos1 = which(out_SKAT_List$param$rho == 1)
        Pvalue = c(out_SKAT_List$p.value,                  # SKAT-O
                   out_SKAT_List$param$p.val.each[pos0],   # SKAT
                   out_SKAT_List$param$p.val.each[pos1])   # Burden Test
        error.code = 0
      }
      
      pval.Region = rbind.data.frame(pval.Region,
                                     data.frame(Region = regionName,
                                                nMarkers = nMarker,
                                                nMarkersURV = nMarkerURV,
                                                Anno.Type = AnnoName,
                                                pval.SKATO = Pvalue[1], 
                                                pval.SKAT = Pvalue[2],
                                                pval.Burden = Pvalue[3]))
    }
    
    # output.Region = cbind.data.frame(info.Region, pval.Region)
    
    # Util.R: write summary statistics to output file.
    # writeOutputFile(output.Region, i, OutputFile, OutputFileIndex, "Region", 1)
    writeOutputFile(list(pval.Region, info.Region), i, 
                    list(OutputFile, paste0(OutputFile, ".markerInfo")), 
                    OutputFileIndex, "Region", 1)
  }
  
  message = paste0("Analysis done! The results have been saved to '", OutputFile,"' and '",
                   paste0(OutputFile, ".markerInfo"),"'.")
  return(message)
}


setRegion = function(NullModelClass, objNull, control, chrom, SparseGRMFile, Group, ifOutGroup)
{
  # The following function is in Main.cpp
  nGroup = length(unique(Group))
  setRegion_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$max_maf_region,
                            control$min_mac_region,
                            control$max_markers_region,
                            control$omp_num_threads,
                            Group, ifOutGroup, nGroup)
  
  # Check POLMM.R
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setRegion = setRegion.POLMM(objNull, control, chrom, SparseGRMFile)  
  
  # To be continued
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setRegion = setRegion.SAIGE(objNull, control, chrom, SparseGRMFile)
  
  # Check SPACox.R
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setRegion = setRegion.SPACox(objNull, control)
  
  return(obj.setRegion)
}

mainRegion = function(NullModelClass, 
                      genoType, 
                      genoIndex, 
                      OutputFile,
                      n,
                      P1Mat,
                      P2Mat)
{
  if(NullModelClass == "POLMM_NULL_Model")
    obj.mainRegion = mainRegionInCPP("POLMM", genoType, genoIndex, OutputFile, n, P1Mat, P2Mat)
  
  if(NullModelClass == "SPACox_NULL_Model")
    obj.mainRegion = mainRegionInCPP("SPACox", genoType, genoIndex, OutputFile, n, P1Mat, P2Mat)
  
  return(obj.mainRegion)
}

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
getRegionList = function(RegionFile,
                         RegionAnnoHeader,
                         markerInfo)
{
  print(paste0("Start extracting marker-level information from 'RegionFile' of ", RegionFile, "."))
  
  if(!file.exists(RegionFile))
    stop(paste("Cannot find 'RegionFile' in", RegionFile))
  
  RegionData = data.table::fread(RegionFile, header = T, stringsAsFactors = F);
  RegionData = as.data.frame(RegionData)
  colnames(RegionData)[1:2] = toupper(colnames(RegionData)[1:2])
  
  if(any(colnames(RegionData)[1:2] != c("REGION", "MARKER")))
    stop("The first two elements in the header of 'RegionFile' should be c('REGION', 'MARKER').")
  
  # updated on 2021-08-05
  colnames(markerInfo)[3] = "MARKER"
  RegionData = merge(RegionData, markerInfo, by = "MARKER", all.x = T, sort = F)
  posNA = which(is.na(RegionData$genoIndex))
  
  if(length(posNA) != 0){
    print(head(RegionData[posNA,1:2]))
    stop(paste0("Total ",length(posNA)," markers in 'RegionFile' are not in 'GenoFile'."))
  }
  
  HeaderInRegionData = colnames(RegionData)
  if(!is.null(RegionAnnoHeader)){
    if(any(!RegionAnnoHeader %in% HeaderInRegionData))
      stop("At least one element in 'RegionAnnoHeader' is not in the header of RegionFile.")
    posAnno = which(HeaderInRegionData %in% RegionAnnoHeader)
  }else{
    print("Since no 'RegionAnnoHeader' is given, region-based testing will not incorporate any annotation information.")
    posAnno = NULL
  }
  
  RegionList = list()
  uRegion = unique(RegionData$REGION)
  for(r in uRegion){
    
    # print(paste0("Analyzing region ",r,"...."))
    
    posSNP = which(RegionData$REGION == r)
    SNP = RegionData$MARKER[posSNP]
    
    if(any(duplicated(SNP)))
      stop(paste0("Please check RegionFile: in region ", r,": duplicated SNPs exist."))
    
    # posMarker = match(SNP, markerInfo$ID, 0)
    # if(any(posMarker == 0))
    #   stop(paste0("At least one marker in region ", r," are not in 'GenoFile' and 'GenoFileIndex'."))
    
    regionMat = cbind(BASE=1, RegionData[posSNP, posAnno, drop=F])
    rownames(regionMat) = SNP
    
    # genoIndex = markerInfo$genoIndex[posMarker]
    # chrom = markerInfo$CHROM[posMarker]
    genoIndex = RegionData$genoIndex[posSNP]
    chrom = RegionData$CHROM[posSNP]
    uchrom = unique(chrom)
    
    if(length(uchrom) != 1)
      stop(paste0("In region ",r,", markers are from multiple chromosomes."))
    
    RegionList[[r]] = list(SNP = SNP,
                           regionMat = regionMat,
                           genoIndex = genoIndex,
                           chrom = uchrom)
  }
  
  return(RegionList)
}