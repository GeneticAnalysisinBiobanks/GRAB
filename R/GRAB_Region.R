
#' GRAB: region-based multiple markers analysis
#' 
#' Region-based analysis: Test for association between multiple variants (mostly low-frequency or rare variants) in a region and an ordinal categorical phenotype via POLMM
#' 
#' @param objNull an output object of the POLMM_Null_Model() function with a class of "POLMM_NULL_Model". 
#' @param AnnoFile a character of annotation file. Column 1: marker ID, Column 2: Region ID, Columns 3-n annotation similar as in STAAR
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". The default is NULL, that is, to share the same prefix as GenoFile. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results
#' @param chrom a character to specify chromosome of the markers in analysis. Must be specified unless LOCO = F when fitting the null model.
#' @param SparseGRM an object of class "SparseGRM", check help(getSparseGRM) for more details.
#' @param POLMM.control a list of parameters for controlling the POLMM.Marker(). The default is NULL, that is, to use the default parameters. Check 'Details' and 'Examples' for more details.
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
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' GenoFile = gsub("-ext.fam", "-ext.bed", famFile)
#' AnnoFile = system.file("extdata", "AnnoFile.txt", package = "POLMM")
#' OutputFile = gsub("AnnoFile","OutputFile",AnnoFile)
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)
#' objNullFile = system.file("objNull.RData", package = "POLMM")
#' load(objNullFile)
#' chrom = 1
#' 
#' OUTPUT = POLMM.Region(objNull, AnnoFile, GenoFile, GenoFileIndex = NULL, OutputFile, 
#'                       SparseGRM, chrom, POLMM.control = list(max_maf_region = 0.5))
#'      
#' @export
#' @import SKAT, data.table

GRAB.Region = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile = NULL,
                       regionFile,              # column 1: marker Set ID, column 2: SNP ID, columns 3-n: Annotations similar as in STAAR
                       regionAnnoHeader = NULL,
                       SparseGRM,
                       chrom,
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);
  checkOutputFile(OutputFile);
  
  ## check the setting of control, if not specified, the default setting will be used
  control = checkControl(control)
  print("The below is the list of control parameters used in analysis.")
  print(control)
  
  SampleInModel = as.character(objNull$subjIDs);
  
  ## add something in setGenoInput(.) to select specific markers requested by users
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleInModel)  
  markers = objGeno$markers
  
  ## annotation
  regionList = getRegionList(regionFile, regionAnnoHeader)
  nRegions = length(regionList)
  
  setRegion(NullModelClass, objNull, control)
  
  for(i in 1:nRegions){
    region = regionList[i]
    regionName = names(region)
    print(paste0("Analyzing Region of ", regionName, "......"))
    
    resRegion = mainRegion(region, NullModelClass, objNull, control)
    
    quote = F
    if(i == 1){
      append = F; col.names = T
    }else{
      append = T; col.names = F
    }
    
    data.table::fwrite(resRegion, OutputFile, quote = quote, sep = "\t", append = append, col.names = col.names)
  }
    
  message = paste0("The analysis results have been saved to '", OutputFile,"'.")
  return(message)
}

setRegion = function()
{
  if(objNull$controlList$LOCO){
    if(!chrom %in% names(objNull$LOCOList))
      stop("'chrom' should be in names(objNull$LOCOList).")
    obj.CHR = objNull$LOCOList[[chrom]]
    
    if(!chrom %in% names(SparseGRM))
      stop("'chrom' should be in names(SparseGRM).")
    SparseGRM.CHR = SparseGRM[[chrom]]
  }
  
  SPmatR.CHR = makeSPmatR(SparseGRM.CHR, objNull$subjIDs)
  
  setPOLMMobjInR(obj.CHR$muMat,
                 obj.CHR$iRMat,
                 objNull$Cova,
                 objNull$yVec,          # 1 to J
                 SPmatR.CHR,
                 objNull$tau,
                 POLMM.control$printPCGInfo,
                 POLMM.control$tolPCG,
                 POLMM.control$maxiterPCG)
  
  memory_chunk = POLMM.control$memory_chunk
  
  # to be continued
  n = length(SubjID.step1)
  J = max(objNull$yVec)
  p = ncol(objNull$Cova)
  NonZero_cutoff = floor(log(1e7, J))  # for efficient resampling (ER)
  maxMarkers = getMaxMarkers(memory_chunk, n, J, p);
  
  print(paste0("The current POLMM.control$memory_chunk is ", memory_chunk,"(GB)."))
  print(paste0("Based on the sample size, we divide region with more than ", maxMarkers, " markers into multiple chunks to save memory usage."))
  print("If the memory usage still exceed the memory you request, please set a smaller POLMM.control$memory_chunk.")
  
  StdStat_cutoff = POLMM.control$SPA_cutoff;
}

# calculate region-based p-values for any given annotations
mainRegion = function(region, NullModelClass, objNull, control)
{
  markers = region$markers
  annoMat = region$annoMat
  
  OutList = mainRegioninCPP(markers, NullModelClass)
  
  ## extract information from mainRegioninCPP(.)
  StatVec = OutList$StatVec
  VarSVec = diag(OutList$VarSMat)
  pvalNormVec = OutList$pvalNormVec;
  adjPVec = OutList$pvalVec;
  adjVarSVec = StatVec^2 / qchisq(adjPVec, df = 1, lower.tail = F)
  weightVec = OutList$weightVec
  
  # index of SNPs passing criterion (MAF, missing rate, et al.)
  posVec = OutList$posVec + 1   # "+1" because C++ starts from 0 and R starts from 1
  annoMat = annoMat[posVec,,drop=F]
  q = ncol(annoMat)             # number of annotations: column 1 is always 1s
  
  r0 = adjVarSVec / VarSVec 
  r0 = pmax(r0, 1)
  
  pvalVec.BT = c()
  pvalVec.SKAT = c()
  pvalVec.SKATO = c()
  annoVec = c()
  errVec = c()
  
  # cycle for multiple annotations
  for(i in 1:q){
    annoName = colnames(annoMat)[i]
    annoWeights = weightVec * annoMat[,i]
    wr0 = sqrt(r0) * annoWeights
    wStatVec = StatVec * annoWeights
    wadjVarSMat = t(OutList$VarSMat * wr0) * wr0
    wadjVarSMat = wadjVarSMat * max(OutList$rBT, 1)
    
    out_SKAT_List = try(SKAT:::Met_SKAT_Get_Pvalue(Score = wStatVec, 
                                                   Phi = wadjVarSMat,
                                                   r.corr = control$r_corr, 
                                                   method = "optimal.adj", 
                                                   Score.Resampling = NULL),
                        silent = TRUE)
    
    # betaVec = StatVec / adjVarSVec; 
    if(class(out_SKAT_List) == "try-error"){
      Pvalue = c(NA, NA, NA)
      errCode = 2
    }else if(!any(c(0,1) %in% out_SKAT_List$param$rho)){
      Pvalue = c(NA, NA, NA)
      errCode = 3
    }else{
      pos0 = which(out_SKAT_List$param$rho == 0)
      pos1 = which(out_SKAT_List$param$rho == 1)
      Pvalue = c(out_SKAT_List$p.value,                  # SKAT-O
                 out_SKAT_List$param$p.val.each[pos0],   # SKAT
                 out_SKAT_List$param$p.val.each[pos1])   # Burden Test
      errCode = 0
    }
    
    pvalVec.SKATO = c(pvalVec.SKATO, Pvalue[1])
    pvalVec.SKAT = c(pvalVec.SKAT, Pvalue[2])
    pvalVec.BT = c(pvalVec.BT, Pvalue[3])
    annoVec = c(annoVec, annoName)
    errVec = c(errVec, errCode)
  }
  
  ###################
  
  OUT.Region = rbind(OUT.Region,
                     c(Region, 
                       length(OutList$markerVec), 
                       paste(annoVec, collapse = ","),
                       paste(pvalVec.SKATO, collapse = ","),
                       paste(pvalVec.SKAT, collapse = ","),
                       paste(pvalVec.BT, collapse = ","),
                       paste(errVec, collapse = ","),
                       paste(OutList$markerVec, collapse = ","),
                       paste(OutList$freqVec, collapse = ","),
                       paste(OutList$flipVec, collapse = ","),
                       paste(StatVec, collapse = ","),
                       paste(adjVarSVec, collapse = ","),
                       paste(adjPVec, collapse = ",")))
  
  tempFiles = list.files(path = dirname(OutputFile),
                         pattern = paste0("^",basename(OutputFile),".*\\.bin$"),
                         full.names = T)
  file.remove(tempFiles)
  
  # out_Multi_Set = rnorm(1)
  colnames(OUT.Region) = c("regionName", "nMarkers", "Annotation", "P.SKAT-O", "P.SKAT", "P.Burden",
                           "error.code", "markerInfo", "markerMAF","markerAlleleFlip",
                           "markerStat","markerVarS","markerPvalue")
  
  return(resRegion)
}
