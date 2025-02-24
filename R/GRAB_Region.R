
#' Conduct region-level genetic association testing
#' 
#' Test for association between phenotype of interest and regions including multiple genetic marker (mostly low-frequency or rare variants).
#' 
#' @param objNull the output object of function \code{\link{GRAB.NullModel}}. 
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the prefix is the same as GenoFile. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param OutputFile a character of output file to save the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point. If the program ends unexpectedly, the end point can help \code{GRAB} package understand where to restart the analysis. If \code{NULL} (default), \code{OutputFileIndex = paste0(OutputFile, ".index")}. 
#' @param GroupFile a character of region file to specify region-marker mapping with annotation information. Each region includes two or three rows. Only alphabet, numbers, and [:,_+-] symbols are supported. Columns are separated by 'tab'. 
#' @param SparseGRMFile a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}.
#' @param SampleFile a character of file to include sample information with header.
#' @param MaxMAFVec a character of multiple max MAF cutoffs (comma separated) to include markers for region-level analysis. Default value is \code{"0.05,0.01,0.005"}.
#' @param annoVec a character of multiple annotation groups (comma separated) to include markers for region-level analysis. Default value is \code{"lof,lof:missense,lof:missense:synonymous"}.
#' @param chrom to be continued 
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
#' GroupFile = system.file("extdata", "example.GroupFile.txt", package = "GRAB")
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' 
#' ## make sure the output files does not exist at first
#' file.remove(OutputFile)
#' file.remove(paste0(OutputFile, ".markerInfo"))
#' file.remove(paste0(OutputFile, ".index"))
#' 
#' GRAB.Region(objNull = objNull,
#'             GenoFile = GenoFile,
#'             GenoFileIndex = NULL,
#'             OutputFile = OutputFile,
#'             OutputFileIndex = NULL,
#'             GroupFile = GroupFile,
#'             SparseGRMFile = SparseGRMFile,
#'             MaxMAFVec = "0.01,0.005")
#'             
#' data.table::fread(OutputFile)
#' data.table::fread(paste0(OutputFile,".markerInfo"))
#' data.table::fread(paste0(OutputFile,".otherMarkerInfo"))
#' data.table::fread(paste0(OutputFile,".index"), sep="\t", header=F)
#' 
#' SampleFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' GRAB.Region(objNull = objNull,
#'             GenoFile = GenoFile,
#'             GenoFileIndex = NULL,
#'             OutputFile = OutputFile,
#'             OutputFileIndex = NULL,
#'             GroupFile = GroupFile,
#'             SparseGRMFile = SparseGRMFile,
#'             SampleFile = SampleFile,
#'             control = list(SampleLabelCol = "OrdinalPheno"))
#'             
#' data.table::fread(OutputFile)
#' data.table::fread(paste0(OutputFile,".markerInfo"))
#' data.table::fread(paste0(OutputFile,".otherMarkerInfo"))
#' data.table::fread(paste0(OutputFile,".index"), sep="\t", header=F)
#' 
#' @export
#' @import SKAT, data.table

GRAB.Region = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,
                       GroupFile,             
                       SparseGRMFile = NULL,
                       SampleFile = NULL,
                       MaxMAFVec = "0.01,0.001,0.0005",
                       annoVec = "lof,lof:missense,lof:missense:synonymous",
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
  control = checkControl.Region(control)
  
  textToParse = paste0("control = checkControl.Region.", method, "(control)")
  eval(parse(text = textToParse))
  
  # print control list 
  print("The below is the list of control parameters used in region-level genetic association analysis.")
  print(control)
  
  MaxMAFVec = MaxMAFVec %>% strsplit(split = ",") %>% unlist() %>% as.numeric()
  if(any(is.na(MaxMAFVec)))
    stop("any(is.na(MaxMAFVec)):\t")
  
  MaxMAF = max(MaxMAFVec)
  if(MaxMAF > 0.05) 
    stop("Maximal value of 'MaxMAFVec' should be <= 0.05.")
  control$max_maf_region = MaxMAF
  
  annoVec = annoVec %>% strsplit(split = ",") %>% unlist()
  annoList = annoVec %>% strsplit(split = ":")
  allAnno = annoList %>% unlist() %>% unique()

  subjData = as.character(objNull$subjData);
  n = length(subjData)
  
  # updated on 2022-05-09: record MAC and MAF in each sample group (not marker group), to be continued
  SampleLabelNumber = rep(1, n)
  SampleLabelLevels = NULL
  if(!is.null(SampleFile))
  {
    SampleInfo = data.table::fread(SampleFile)
    if(colnames(SampleInfo)[1] != "IID")
      stop("The header of the first column in 'SampleFile' should be 'IID'.")
    
    pos = which(!subjData %in% SampleInfo$IID)
    if(length(pos) > 0)
      stop("At least one subject in null model fitting does not in 'SampleFile'.\n",
           paste0(subjData[pos], collapse = "\t"))
    
    if(!is.null(control$SampleLabelCol))
    {
      SampleLabelColName = control$SampleLabelCol
      if(!SampleLabelColName %in% colnames(SampleInfo))
        stop("'SampleFile' should include one column with header of:\t", SampleLabelColName)
      
      posInSampleInfo = match(subjData, SampleInfo$IID)
      SampleLabel = SampleInfo[[SampleLabelColName]][posInSampleInfo]
      SampleLabelFactor = as.factor(SampleLabel)
      SampleLabelNumber = as.numeric(SampleLabelFactor)
      SampleLabelLevels = levels(SampleLabelFactor)
    }
  }
  nLabel = max(SampleLabelNumber)
  
  ## set up an object for genotype data
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData, control)  # Check 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo

  RegionList = getInfoGroupFile(GroupFile)
  nRegions = length(RegionList)
  
  with(control, 
       setRegion_GlobalVarsInCPP(impute_method, 
                                 missing_cutoff, 
                                 max_maf_region, 
                                 min_mac_region, 
                                 max_markers_region, 
                                 omp_num_threads,
                                 weights.beta,
                                 MaxMAFVec))
  
  textToParse = paste0("obj.setRegion = setRegion.", method, "(objNull, control, chrom, SparseGRMFile)")
  eval(parse(text = textToParse))
  
  diffTime1 = 0
  diffTime2 = 0
  diffTime3 = 0
  
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
    
    # if(nMarkers == 0)
    #   stop("nrow(regionInfo) == 0: no markers are found for region '", regionID, "'.")
    if(nMarkers == 0)
      next;
    
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
        stop("Marker weights cannot be non-positive (<= 0) or NA.")
    }
    
    print(paste0("Analyzing Region of ", regionID, " (",i,"/",nRegions,")."))
    cat("Total", length(regionInfo$ID), "markers:\t", 
        paste0(head(regionInfo$ID), collapse = ", "), "\n")

    t11 = Sys.time()
    obj.mainRegionInCPP = mainRegionInCPP(method, genoType, genoIndex, weightVec, OutputFile, 
                                          SampleLabelNumber, nLabel, 
                                          annoMat, annoVec)
    t12 = Sys.time()
    diffTime1 = diffTime1 + (t12-t11)
    
    # updated on 2022-06-24 (save sum of genotype to conduct burden test and adjust p-values using SPA)
    pvalBurden = obj.mainRegionInCPP$pvalBurden
    # print(class(pvalBurden))
    
    # updated on 2023-02-06 (record summary statistics for sum of genotype for a region)
    infoBurdenNoWeight = obj.mainRegionInCPP$infoBurdenNoWeight
    infoBurdenNoWeight = as.data.frame(infoBurdenNoWeight)
    infoBurdenNoWeight = cbind(regionID, infoBurdenNoWeight)
    colnames(infoBurdenNoWeight) = c("region", "anno", "max_maf", "sum", "Stat", "beta", "se.beta", "pvalue")
    
    infoBurdenNoWeight$anno = annoVec[infoBurdenNoWeight$anno + 1]
    infoBurdenNoWeight$max_maf = MaxMAFVec[infoBurdenNoWeight$max_maf + 1]
    
    ## add annotation information
    obj.mainRegionInCPP$AnnoVec = c(regionInfo$Annos, annoVec)
    if(!is.null(SampleLabelLevels))
    {
      colnames(obj.mainRegionInCPP$MACLabelMat) = paste0("MAC_", SampleLabelLevels)
      colnames(obj.mainRegionInCPP$MAFLabelMat) = paste0("MAF_", SampleLabelLevels)
    }
    
    textToParse = paste0("obj.mainRegion = mainRegion.", method, "(genoType, genoIndex, OutputFile, control, n, obj.setRegion, obj.mainRegionInCPP, nLabel)")
    eval(parse(text = textToParse))
    
    Other.Markers = obj.mainRegion$Other.Markers %>% mutate(Region = regionID, .before = ID)
    VarMat = obj.mainRegion$VarMat
    RV.Markers0 = obj.mainRegion$RV.Markers %>% mutate(Region = regionID, .before = ID)
    
    Other.Markers = regionInfo %>% select(ID, Annos) %>% merge(Other.Markers, by = "ID")
    
    if(nrow(VarMat) != nrow(RV.Markers0))
      stop("nrow(VarMat) != nrow(RV.Markers0)!")
    
    RV.Markers = RV.Markers0 %>% 
      mutate(betaWeights = dbeta(MAF, control$weights.beta[1], control$weights.beta[2]),
             adjVarSVec = StatVec^2 / qchisq(pval1Vec, df = 1, lower.tail = F),
             # r0 = adjVarSVec / diag(VarMat),  # edited on 06/22/2022
             r0 = pmax(adjVarSVec / diag(VarMat), 1),
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
    
    t21 = Sys.time()
    pval.Region = data.frame()
    iSPA = 1
    for(anno in annoVec)
    {
      annoTemp = unlist(strsplit(anno, split = ":"))
      
      posURV = RV.MarkersURV %>% filter(ID == anno) %>% select(posRow) %>% unlist()
      nMarkersURV = Other.MarkersWithAnno %>% filter(Annos %in% annoTemp) %>% nrow()
      if(length(posURV) != 1)
        stop("length(posURV) != 1")
      
      for(MaxMAF in MaxMAFVec)
      {
        posRV = RV.MarkersWithAnno %>% filter(MAF < MaxMAF & Annos %in% annoTemp) %>% select(posRow) %>% unlist()
        pos = c(posRV, posURV)
        n1 = length(pos)
        
        ScoreBurden = sum(RV.Markers$wStatVec[pos])
        VarBurden = sum(wadjVarSMat[pos, pos])
        pvalBurdenSPA = pvalBurden[iSPA,2]
        VarBurdenSPA = ScoreBurden^2 / qchisq(pvalBurdenSPA, df = 1, lower.tail = F)
        ratioBurdenSPA = max(VarBurdenSPA/VarBurden, 1)
        iSPA = iSPA + 1
        
        t31 = Sys.time()
        out_SKAT_List = with(RV.Markers, try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec[pos], 
                                                                        # Phi = wadjVarSMat[pos, pos],  
                                                                        Phi = ratioBurdenSPA * wadjVarSMat[pos, pos],  
                                                                        r.corr = control$r.corr, 
                                                                        method = "optimal.adj", 
                                                                        Score.Resampling = NULL),
                                             silent = TRUE))

        t32 = Sys.time()
        diffTime3 = diffTime3 + (t32-t31)
        
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
    
    ## Cauchy Combination
    pval.Cauchy.SKATO = CCT(pval.Region$pval.SKATO)
    pval.Cauchy.SKAT = CCT(pval.Region$pval.SKAT)
    pval.Cauchy.Burden = CCT(pval.Region$pval.Burden)
    
    pval.Region = rbind.data.frame(pval.Region,
                                   data.frame(Region = regionID,
                                              nMarkers = NA,
                                              nMarkersURV = NA,
                                              Anno.Type = "Cauchy",
                                              MaxMAF.Cutoff = NA,
                                              pval.SKATO = pval.Cauchy.SKATO, 
                                              pval.SKAT = pval.Cauchy.SKAT,
                                              pval.Burden = pval.Cauchy.Burden))
    
    t22 = Sys.time()
    diffTime2 = diffTime2 + (t22-t21)
    
    # cat("diffTime1:\t", diffTime1,"\n")
    # cat("diffTime2:\t", diffTime2,"\n")
    # cat("diffTime3:\t", diffTime3,"\n")
    
    writeOutputFile(Output = list(pval.Region, 
                                  RV.Markers0, 
                                  Other.Markers,
                                  infoBurdenNoWeight), 
                    OutputFile = list(OutputFile, 
                                      paste0(OutputFile, ".markerInfo"), 
                                      paste0(OutputFile, ".otherMarkerInfo"),
                                      paste0(OutputFile, ".infoBurdenNoWeight")), 
                    OutputFileIndex = OutputFileIndex,
                    AnalysisType = "Region",
                    nEachChunk = 1,
                    indexChunk = i,
                    Start = (i==1),
                    End = (i==nRegions))
  }
      
  cat("mainRegionInCPP():\t", diffTime1,"seconds.\n")
  # cat("SKATO:\t", diffTime2,"seconds.\n")
  cat("SKATO::\t", diffTime3,"\n")
  printTimeDiffInCPP()
  
  cat("Analysis done! The results have been saved to the below files:\n", 
      OutputFile, "\n",
      paste0(OutputFile, ".markerInfo\n"),
      paste0(OutputFile, ".otherMarkerInfo\n"))
  
  message = paste0("Analysis done! The main results have been saved to '", OutputFile,"'")
  
  return(message)
}

 
getInfoGroupLine = function(markerGroupLine, nLine)
{
  if(length(markerGroupLine) == 0)
    stop("The line ", nLine," in `groupFile` is empty.")
  
  info = markerGroupLine %>% strsplit(split="\t") %>% unlist()
  if(length(info) < 3)
    stop("The line ", nLine, " in 'groupFile' includes < 3 elements, please note that each line should be seperated by 'tab'.")
  
  geneID = info[1];
  type = info[2];
  values = info[c(-1,-2)]
  
  grepTemp = grep(" ", values, value = T)
  if(length(grepTemp) > 0)
    stop("'GroupFile' cannot contain 'space':\n", 
         grepTemp %>% unique() %>% paste0(collapse = "\t"))
  
  grepTemp = grep(";", values, value = T)
  if(length(grepTemp) > 0)
    stop("'GroupFile' cannot contain ';':\n", 
         grepTemp %>% unique() %>% paste0(collapse = "\t"))

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
      if(nRegion == 1)
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
  
    if(!type %in% c("var", "anno", "weight"))
      stop("The second column of the groupFile (tab-seperated) should be one of 'var', 'anno', and 'weight'.\n
         Please double check line ", nLine, ".")
    
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
