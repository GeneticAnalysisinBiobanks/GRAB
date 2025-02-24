#' Quality control to check batch effect between study cohort and reference population.
#'
#' Quality control to check batch effect between study cohort and reference population.
#' 
#' @param GenoFile a character of genotype file. See Details section for more details.
#' @param GenoFileIndex additional index file(s) corresponding to GenoFile. See Details section for more details.
#' @param PhenoData a dataframe. must have at least two columns \code{SampleID} and \code{Indicator}, \code{SampleID} hold the personal identifiers for all individuals, and \code{Indicator} hold wether the event occurred, and the value must be 0 or 1 or NA 
#' @param RefAFfile a character of reference file. The reference file must be a \code{txt} file (header are required) including at least 7 columns \code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}, \code{AF_ref}, \code{AN_ref}  
#' @param RefPrev a numeric of the event rate in the population
#' @param SNPnum a int. The least number of markers. The default is 1e4
#' @param control a list of parameters to decide which markers to extract.See \code{Details} section for more details
#' @return an R object with a class of "QCforBatchEffect".
#' \itemize{
#'   \item{mergeGenoInfo}: a dataframe of marker info and reference MAF
#'   \item{cutoff}: a numeric, the cut off of batcheffect.
#'   \item{count}: a dataframe of the frequency of the batch effect pvalue.
#'   \item{PhenoData}: a dataframe of the input PhenoData.
#'   \item{control}: a list of parameters to decide which markers to extract.
#' }
#' @details
#' ## The following details are about argument \code{control}
#' \describe{
#' Argument \code{control} is used to include and exclude markers. The function supports two include files of (\code{IDsToIncludeFile}, \code{RangesToIncludeFile}) and two exclude files of (\code{IDsToExcludeFile}, \code{RangesToExcludeFile}), but does not support both include and exclude files at the same time. 
#'   \itemize{
#'     \item \code{AlleleOrder}: a character, \code{"ref-first"} or \code{"alt-first"}, to determine whether the REF/major allele should appear first or second. 
#'     \item \code{AllMarkers}:  a logical value (default: \code{FALSE}) to indicate if all markers are extracted. It might take too much memory to put genotype of all markers in R. 
#'     \item \code{IDsToIncludeFile}: a file of marker IDs to include, one column (no header). Check \code{system.file("extdata", "IDsToInclude.txt", package = "GRAB")} for an example.
#'     \item \code{IDsToExcludeFile}: a file of marker IDs to exclude, one column (no header).
#'     \item \code{RangesToIncludeFile}: a file of ranges to include, three columns (no headers): chromosome, start position, end position. Check \code{system.file("extdata", "RangesToInclude.txt", package = "GRAB")} for an example.
#'     \item \code{RangesToExcludeFile}: a file of ranges to exclude, three columns (no headers): chromosome, start position, end position.
#'   }
#' }
#' 
#' 
#' @export
#' @import dplyr, data.table
#' @examples
#' Check ?GRAB.WtSPAG for examples
QCforBatchEffect = function(GenoFile               # a character of file names of genotype files
                            ,GenoFileIndex = NULL  # additional index file(s) corresponding to GenoFile
                            ,OutputFile
                            ,control=list(AlleleOrder = "ref-first")  
                            ,PhenoData             # an R data frame with at least two columns, headers are required and should include c("SampleID", "Indicator"), the "Indicator" column should be 0, 1, or NA. 
                            ,RefAfFile             # a character of file name of refInfo, which including at least 7 columns, header are required and should include c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref","AN_ref")  
                            ,RefPrevalence         # refernce population prevalence, the proportion of indicator == 1.
                            ,SNPnum=1e4            # default least number of SNPs is 1e4
                            
){
  if(is.null(OutputFile))
    stop("Argument of 'OutputFile' is required to store information for the follow-up analysis.")
  
  # check if there are c("Indicator", "SampleID") in PhenoData-------------------
  if(!is.null(control$IndicatorColumn))
  {
    if(!control$IndicatorColumn %in% colnames(PhenoData))
      stop(paste0("Cannot find a column of '",
                  control$IndicatorColumn,
                  "' (i.e. control$IndicatorColumn) in colnames(PhenoData)"))
    posCol = which(colnames(PhenoData) == control$IndicatorColumn)
    colnames(PhenoData)[posCol] = "Indicator"
  }
  
  if(!is.null(control$SampleIDColumn))
  {
    if(!control$SampleIDColumn %in% colnames(PhenoData))
      stop(paste0("Cannot find a column of '",
                  control$IndicatorColumn,
                  "' (i.e. control$SampleIDColumn) in colnames(PhenoData)"))
    posCol = which(colnames(PhenoData) == control$SampleIDColumn)
    colnames(PhenoData)[posCol] = "SampleID"
  }  
  
  if(!"Indicator" %in% colnames(PhenoData))
    stop("The column of 'Indicator' is required in PhenoData!")
  
  if(any(!unique(PhenoData$Indicator) %in% c(0,1,NA)))
    stop("The value of Indicator should be 0,1 or NA")
  
  if(!"SampleID" %in% colnames(PhenoData))
    stop("The column of 'SampleID' is required in PhenoData!")
  
  #step1: quality control--------------------------------------------------------
  suppressPackageStartupMessages(library("GRAB",quietly = T))
  suppressPackageStartupMessages(library("data.table",quietly = T))
  suppressPackageStartupMessages(library("dplyr",quietly = T))
  ## reference genoInfo----------------------------------------------------------
  refGenoInfo = fread(RefAfFile)%>%as_tibble()
  
  # check if there are 7 columns in RefAfFile
  for(colname in c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref","AN_ref")){
    if(!colname %in% colnames(refGenoInfo)){
      stop( paste0(colname, " is missing in RefAfFile!") )}
  }
  # ## merge sample genoInfo and ref genoInfo--------------------------------------
  # GenoInfo.ctrl = GRAB.getGenoInfo(GenoFile = GenoFile
  #                                  ,GenoFileIndex = GenoFileIndex
  #                                  ,SampleIDs = with(PhenoData,SampleID[Indicator==0])
  #                                  ,control = control)%>%
  #   rename(mu0 = altFreq, mr0 = missingRate ) %>% select(mu0, mr0)
  # mergeGenoInfo = GRAB.getGenoInfo(GenoFile = GenoFile
  #                                  ,GenoFileIndex = GenoFileIndex
  #                                  ,SampleIDs = with(PhenoData,SampleID[Indicator==1])
  #                                  ,control = control) %>% 
  #   rename(mu1 = altFreq , mr1 = missingRate) %>% 
  #   cbind(., GenoInfo.ctrl)%>%as_tibble()%>%
  #   merge(.,refGenoInfo,by=c("CHROM", "POS", "REF", "ALT", "ID"), all.x=T)
  
  ## merge sample genoInfo and ref genoInfo--------------------------------------
  GenoInfo.ctrl = GRAB.getGenoInfo(GenoFile = GenoFile
                                   ,GenoFileIndex = GenoFileIndex
                                   ,SampleIDs = with(PhenoData,SampleID[Indicator==0])
                                   ,control = control)%>%
    rename(mu0 = altFreq, mr0 = missingRate ) %>% select(mu0, mr0)
  
  
  GenoInfo = GRAB.getGenoInfo(GenoFile = GenoFile
                              ,GenoFileIndex = GenoFileIndex
                              ,SampleIDs = with(PhenoData,SampleID[Indicator==1])
                              ,control = control) %>% 
    rename(mu1 = altFreq , mr1 = missingRate) %>% 
    cbind(., GenoInfo.ctrl)%>%as_tibble()%>%
    mutate(RA = ifelse(REF<ALT, paste0(REF, ALT), paste0(ALT,REF)))
  
  mergeGenoInfo = refGenoInfo %>% 
    mutate(RA = ifelse(REF<ALT, paste0(REF, ALT), paste0(ALT,REF))) %>% 
    merge(.,GenoInfo,by=c("CHROM", "POS", "RA"),all.y=T) %>% 
    rename(REF = REF.y, ALT= ALT.y, ID = ID.y)%>% 
    mutate(AF_ref = ifelse(REF == REF.x, AF_ref , 1-AF_ref  ))%>%
    select(-REF.x,-ALT.x, -ID.x, -RA)
  
  
  ## evaluate batch effect and calculate pvalue----------------------------------
  test = function(n0,n1,n2,w0,w1, p1,p2){
    p = sum(p1 *(n1* w1+n0* w0) + p2*n2)/sum(n1* w1+n0* w0+ n2)
    v=((n1*w1^2+n0*w0^2)/(2*(n1*w1+n0*w0)^2) + 1/(2*n2))*p*(1-p)
    z = (p1 - p2) / sqrt(v)
    p= 2*pnorm(-abs(z), lower.tail=TRUE)
    chisq = z^2
    return(p)
  }
  
  n1=sum(PhenoData$Indicator)*(1-mergeGenoInfo$mr1)
  n0=sum(1-PhenoData$Indicator)*(1-mergeGenoInfo$mr0)
  w1=1; w0=(1-RefPrevalence)/RefPrevalence*n1/n0
  
  pvalue_bat = lapply(1: nrow(mergeGenoInfo), function(ind){
    w.maf = with(mergeGenoInfo, sum(mu0[ind]*w0[ind]*n0[ind] +mu1[ind]*w1*n1[ind] )/(n0[ind]*w0[ind]+n1[ind]*w1))
    p.test = test(n0=n0[ind],n1=n1[ind],n2 =mergeGenoInfo$AN_ref[ind]/2, w0=w0[ind], w1=w1, p1= w.maf, p2=mergeGenoInfo$AF_ref[ind])
  })%>%unlist()
  mergeGenoInfo = mergeGenoInfo %>%mutate(pvalue_bat)
  
  count = table(cut(na.omit(pvalue_bat), breaks = seq(0,1,0.01)))%>%
    as.data.frame() %>% setNames(c("interval", "Freq")) %>%
    mutate(breaks =  head(seq(0,1,0.01), -1))
  
  ## CUT OFF---------------------------------------------------------------------
  if(nrow(mergeGenoInfo)<SNPnum){
    warning(paste0("Please input at least ", SNPnum, " SNPs"))
    cutoff = NA
  }else{
    cutoff = getCutoff(count)
  }
  
  data.table::fwrite(mergeGenoInfo, OutputFile, row.names = F, col.names = T, quote = F, sep = "\t")
  
  return(list(# mergeGenoInfo=mergeGenoInfo, 
    cutoff=cutoff, 
    count=count,
    PhenoData = PhenoData ,
    control = control))
}

getCutoff = function(count){
  right_mean= lapply(2:nrow(count)-1,function(i){
    m = count[(i+1):nrow(count), "Freq"]%>%mean()
    return(m)
  }) %>% unlist()
  
  for(i in 2:nrow(count)-1){
    diff = (count$Freq[i]-right_mean[i])/right_mean[i]
    cat("The diff between ",i,"th interval with the rest intervals:"
        , diff,"\n")
    if(abs(diff)<0.1){
      cutoff = count$breaks[i]
      cat("cutoff=", cutoff,"\n")
      break
    }
  }
  return(cutoff)
}



#' Merge quality control results from multiple genotype files.
#'
#' Merge quality control results from multiple genotype files.
#' 
#' @param qcResult1 the output R list of function QCforBatchEffect().
#' @param ... other output R lists of function QCforBatchEffect().
#' @return an R object with a class of "QCforBatchEffect".
#' \itemize{
#'   \item{cutoff}: a numeric. A new cutoff calculated by using merged QcResults .
#'   \item{allcount}: a dataframe. Record the frequency of all the batch effect pvalue from \code{qcResults}
#' }
#' 
mergeQCresults = function(qcResult1,...){ #  allQcResults is a list including all the QC results from QCforBatchEffect
  
  allQcResults = list(qcResult1,...)
  #check wehter the inputs of QCforBatchEffect function are consistent
  QcInfo = lapply(allQcResults, function(obj){
    n1 = sum(obj$PhenoData$Indicator)
    n0 = sum(1 - obj$PhenoData$Indicator)
    AllelOrder = obj$control$AlleleOrder
    return(cbind(n1,n0,AllelOrder))
  })%>%do.call("rbind" ,.) %>%as_tibble()%>%setNames(c("n1","n0","AllelOrder"))
  if(length(unique(QcInfo$n1))!=1 | length(unique(QcInfo$n0))!=1){
    stop("PhenoData in QCforBatchEffect should be consistent!")
  }
  
  if(length(unique(QcInfo$AllelOrder ))!=1 ){
    stop("AllelOrder in QCforBatchEffect should be consistent!")
  }
  
  
  #ref = lapply(allQcResults, function(q){q$control})
  # Merge all counts
  allFreq=rep(0, length(allQcResults[[1]]$count$Freq))
  for(i in 1: length(allQcResults)){
    allFreq = allFreq + allQcResults[[i]]$count$Freq
  }
  allcount = allQcResults[[1]]$count %>% mutate(Freq = allFreq)
  
  #calculate new cutoff 
  cutoff=getCutoff(allcount)
  return(list(cutoff=cutoff, allcount=allcount ))
}
