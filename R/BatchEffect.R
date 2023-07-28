#' Quality control to check batch effect between study cohort and reference population.
#'
#' Quality control to check batch effect between study cohort and reference population.
#' 
#' @param GenoFile XXX.
#' @param GenoFileIndex XXX.
#' @param PhenoData XXX.
#' @param RefAFfile XXX.
#' @param RefPrev XXX.
#' @param SNPnum XXX.
#' @param control XXX.
#' @return an R object with a class of "QCforBatchEffect".
#' \itemize{
#'   \item{mergeGenoInfo}: XXX.
#'   \item{cutoff}: XXX.
#'   \item{count}: XXX.
#'   \item{PhenoData}: XXX.
#'   \item{control}: XXX.
#' }
#' @details
#' something
#'  
#' \describe{
#' something
#' }
#' 
#' ## The following details are about argument \code{control}
#' \describe{
#' Argument \code{control} includes a list of parameters for controlling the null model fitting process.  
#'   \itemize{
#'     \item \code{maxiter}: Maximum number of iterations used to fit the null model. *(default=100)*
#'     \item \code{seed}: An integer as a random seed. Used when random process is involved. *(default=12345678)*
#'     \item \code{tolBeta}: Positive tolerance: the iterations converge when |beta - beta_old| / (|beta| + |beta_old| + tolBeta) < tolBeta. *(default=0.001)*
#'     \item \code{showInfo}: Whether to show more detailed information for trouble shooting. *(default=FALSE)*
#'   }
#' }
#' 
#' @examples
#' # something of comments followed by R codes
#' load("/gdata01/user/liying/UKB/AD_step2/chr21_1_step2.RData")
#' PhenoData = PheData1[,c("ID","PC1","PC2","PC3","time","event")]%>%rename( SampleID =ID,Indicator=event )
#' 
#' qcResult1 = QCforBatchEffect(GenoFile = c("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c21_b0_v3.bgen")
#'                              ,GenoFileIndex = c("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr21_v3.bgen.bgi", 
#'                                                 "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample")
#'                              ,control=list(AlleleOrder = "ref-first", IDsToIncludeFile = list.files("/gdata02/master_data1/UK_Biobank/ukb22828_imp/split_R2_0.6_MAF_5e-04",paste0("chr_",21,"_chunk_",1,".txt"),full.names = T))
#'                              ,PhenoData=PhenoData
#'                              ,RefAfFile = "/gdata01/user/liying/gnomAD/uk10k.txt"
#'                              ,RefPrevalence = 0.02
#'                              ,SNPnum=1e4
#' )
#' 
#' 
#' qcResult2 = QCforBatchEffect(GenoFile = c("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c22_b0_v3.bgen")
#'                              ,GenoFileIndex = c("/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_imp_chr22_v3.bgen.bgi", 
#'                                                 "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample")
#'                              ,control=list(AlleleOrder ="ref-first", IDsToIncludeFile = list.files("/gdata02/master_data1/UK_Biobank/ukb22828_imp/split_R2_0.6_MAF_5e-04",paste0("chr_",22,"_chunk_",1,".txt"),full.names = T))
#'                              ,PhenoData=PhenoData
#'                              ,RefAfFile = "/gdata01/user/liying/gnomAD/uk10k.txt"
#'                              ,RefPrevalence = 0.02
#'                              ,SNPnum=1e4
#' )
#' 
#' allQC = mergeQCresults(qcResult1,qcResult2)
#' Residual = getResid.wCox(formula = Surv(time, Indicator) ~ .
#'                          , PhenoData=PhenoData
#'                          , RefPrevalence=0.02)
#' @export
#' @import dplyr, data.table
QCforBatchEffect = function(GenoFile              #a character of file names of genotype files
                            ,GenoFileIndex         #additional index file(s) corresponding to GenoFile
                            ,control=list(AlleleOrder = "ref-first")  
                            ,PhenoData             #a dataframe with at least two columns, headers are required and should include c("SampleID", "Indicator"), the "Indicator" column should be 0, 1, or NA. 
                            ,RefAfFile             #a character of file name of refInfo, which including at least 7 columns, header are required and should include c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref","AN_ref")  
                            ,RefPrevalence         #refernce population prevalence, the proportion of indicator == 1.
                            ,SNPnum=1e4            #default least number of SNPs is 1e4
                            
){
  # check if there are c("Indicator", "SampleID") in PhenoData-------------------
  if(!"Indicator" %in% colnames(PhenoData))
    stop("Indicator is missing in PhenoData!")
  
  if(any(!unique(PhenoData$Indicator) %in% c(0,1,NA)))
    stop("The value of Indicator should be 0,1 or NA")
  
  if(!"SampleID" %in% colnames(PhenoData))
    stop("SampleID is missing in PhenoData!")
  
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
  ## merge sample genoInfo and ref genoInfo--------------------------------------
  GenoInfo.ctrl = GRAB.getGenoInfo(GenoFile = GenoFile
                                   ,GenoFileIndex = GenoFileIndex
                                   ,SampleIDs = with(PhenoData,SampleID[Indicator==0])
                                   ,control = control)%>%
    rename(mu0 = altFreq, mr0 = missingRate ) %>% select(mu0, mr0)
  mergeGenoInfo = GRAB.getGenoInfo(GenoFile = GenoFile
                                   ,GenoFileIndex = GenoFileIndex
                                   ,SampleIDs = with(PhenoData,SampleID[Indicator==1])
                                   ,control = control) %>% 
    rename(mu1 = altFreq , mr1 = missingRate) %>% 
    cbind(., GenoInfo.ctrl)%>%as_tibble()%>%
    merge(.,refGenoInfo,by=c("CHROM", "POS", "REF", "ALT","ID"),all.x=T)
  
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
  
  return(list(mergeGenoInfo=mergeGenoInfo, 
              cutoff=cutoff, 
              count=count,
              PhenoData = PhenoData ,
              control = control))
}

#' Merge quality control results from multiple genotype files.
#'
#' Merge quality control results from multiple genotype files.
#' 
#' @param qcResult1 the output R list of function QCforBatchEffect().
#' @param ... other output R lists of function QCforBatchEffect().
#' @return an R object with a class of "QCforBatchEffect".
#' \itemize{
#'   \item{mergeGenoInfo}: XXX.
#'   \item{cutoff}: XXX.
#'   \item{count}: XXX.
#'   \item{PhenoData}: XXX.
#'   \item{control}: XXX.
#' }
mergeQCresults = function(qcResult1,...){ #  allQcResults is a list including all the QC results from QCforBatchEffect
  
  allQcResults = list(...)
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

#' Get model residuals using a weighted Cox regression model.
#'
#' Get model residuals using a weighted Cox regression model.
#' 
#' @param formula an R formula
#' @param PhenoData a dataframe including at least 3 columns c("SampleID","Indicator","time")
#' @param RefPrevalence consisitent with the RefPrevalence used in QCforBatchEffect.
#' @return an R object with a class of "QCforBatchEffect".
#' \itemize{
#'   \item{mergeGenoInfo}: XXX.
#'   \item{cutoff}: XXX.
#'   \item{count}: XXX.
#'   \item{PhenoData}: XXX.
#'   \item{control}: XXX.
#' }
getResid.wCox = function(formula
                         ,PhenoData # a dataframe including at least 3 columns c("SampleID","Indicator","time")
                         ,RefPrevalence # consisitent with the RefPrevalence used in QCforBatchEffect
){
  suppressPackageStartupMessages(library("survival",quietly = T))
  
  # check 
  if(!"Indicator" %in% colnames(PhenoData))
    stop("Indicator is missing in PhenoData!")
  
  if(any(!unique(PhenoData$Indicator) %in% c(0,1,NA)))
    stop("The value of Indicator should be 0,1 or NA")
  
  if(!"SampleID" %in% colnames(PhenoData))
    stop("SampleID is missing in PhenoData!")
  
  if(!"time" %in% colnames(PhenoData))
    stop("time is missing in PhenoData!")
  
  #calculate weight
  PhenoData = PhenoData %>% mutate(weight = ifelse(PhenoData$Indicator == 1,1,(1-RefPrevalence)/RefPrevalence*sum(PhenoData$Indicator)/sum(1-PhenoData$Indicator)))
  #fit null model
  obj.null = coxph(formula, PhenoData, weight = weight, robust = T)
  Residual = data.frame(Resid = obj.null$residuals,SampleID = PhenoData$SampleID )
  return(Residual)
}
