
#' GRAB: region-based multiple markers analysis
#' 
#' Region-based analysis: Test for association between multiple variants (mostly low-frequency or rare variants) in a region and an ordinal categorical phenotype via POLMM
#' 
#' @param objNull an output object of the POLMM_Null_Model() function with a class of "POLMM_NULL_Model". 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". The default is NULL, that is, to share the same prefix as GenoFile. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results
#' @param RegionFile a character of annotation file. Column 1: region ID, Column 2: marker ID, Columns 3-n annotation (non-negative) similar as in STAAR. The header is required and the first two should be "REGION" and "MARKER".
#' @param RegionAnnoHeader a character vector of annotation in analysis. Optional, if not specified, all annotation columns will be used. 
#' @param control a character to specify chromosome of the markers in analysis. Must be specified unless LOCO = F when fitting the null model.
#' @return an R matrix with the following elements
#' \item{ID}{Marker IDs from colnames(GMat)}
#' \item{chr}{Chromosome name from chrVec}
#' \item{MAF}{MAFs of the markers}
#' \item{missing.rate}{Missing rates of the markers}
#' \item{Stat}{Score statistics}
#' \item{VarW}{Estimated variance (VarW) from non-mixed model}
#' \item{VarP}{Estimated variance after adjusting for variance ratio r (VarP = VarW * r)}
#' \item{beta}{Estimated effect size: Stat / VarP}
#' \item{pval.norm}{p values calculated from normal approximation}
#' \item{pval.spa}{p values calculated from saddlepoint approximation}
#' \item{switch.allele}{a logical value indicating if the REF/ALT alleles were switched, if AF > 0.5, we use GVec = 2-GVec, and then give switch.allele=T. This is useful to estimate the effect direction.}
#' @details 
#' More information about the list of 'SKAT.control'
#' \itemize{
#' \item{memory_chunk: a cutoff (Gb) to determine how many markers are in one chunk for region-based analysis [default=4].}
#' \item{kernel: how to weight markers for region-based analysis [default="linear.weighted"].}
#' \item{method: method to conduct region-based analysis [default="method"].}
#' \item{weights_beta:  [default=c(1,25)].}
#' \item{weights:  [default=NULL].}
#' \item{r_corr:  [default=NULL].}
#' \item{printPCGInfo:  [default=FALSE].}
#' \item{tolPCG:  [default=1e-5].}
#' \item{maxiterPCG:  [default=100].}
#' }
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
    
  for(i in outIndex:nRegions){
    
    region = RegionList[[i]]
    regionName = names(RegionList)[i]
    
    SNP = region$SNP
    regionMat = region$regionMat  # annotation values
    genoIndex = region$genoIndex
    chrom = region$chrom
    
    print(paste0("Analyzing Region of ", regionName, "......"))
    print(paste(SNP, collapse = ", "))
    
    obj.setRegion = setRegion(NullModelClass, objNull, control, chrom, SparseGRMFile)
    
    # main function to calculate summary statistics for region-based analysis 
    obj.mainRegion = mainRegion(NullModelClass, genoType, genoIndex, OutputFile, n)
    
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
    
    info.Region = data.frame(Region = regionName,
                             Marker = obj.mainRegion$markerVec,
                             Info = obj.mainRegion$infoVec,
                             AltFreq = obj.mainRegion$altFreqVec,
                             MissingRate = obj.mainRegion$missingRateVec,
                             Beta = obj.mainRegion$BetaVec,
                             seBeta = obj.mainRegion$seBetaVec,
                             pval0 = obj.mainRegion$pval0Vec,
                             pval1 = obj.mainRegion$pval1Vec)
    
    pval.Region = data.frame()
    
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
                                     data.frame(Anno.Type = AnnoName,
                                                pval.SKATO = Pvalue[1], 
                                                pval.SKAT = Pvalue[2],
                                                pval.Burden = Pvalue[3]))
    }
    
    # output.Region = cbind.data.frame(info.Region, pval.Region)
    
    # Util.R: write summary statistics to output file.
    # writeOutputFile(output.Region, i, OutputFile, OutputFileIndex, "Region", 1)
    writeOutputFile(list(pval.Region, info.Region), i, 
                    c(OutputFile, paste0(OutputFile, ".markerInfo")), 
                    OutputFileIndex, "Region", 1)
  }
  
  message = paste0("The analysis results have been saved to '", OutputFile,"' and '",
                   paste0(OutputFile, ".markerInfo"),"'.")
  return(message)
}


setRegion = function(NullModelClass, objNull, control, chrom, SparseGRMFile)
{
  # The following function is in Main.cpp
  setRegion_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$max_maf_region,
                            control$min_mac_region,
                            control$max_mem_region,
                            control$omp_num_threads)
  
  # The following function is in POLMM.R
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setRegion = setRegion.POLMM(objNull, control, chrom, SparseGRMFile)  # check POLMM.R
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setRegion = setRegion.SAIGE(objNull, control, chrom, SparseGRMFile)
  
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setRegion = setRegion.SPACox(objNull, control, chrom)
  
  return(obj.setRegion)
}

mainRegion = function(NullModelClass, 
                      genoType, 
                      genoIndex, 
                      OutputFile,
                      n)
{
  if(NullModelClass == "POLMM_NULL_Model")
    obj.mainRegion = mainRegionInCPP("POLMM", genoType, genoIndex, OutputFile, n)
  
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