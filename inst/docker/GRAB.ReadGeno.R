
# cd /home/wenjianb/Docker/GRAB
# Rscript GRAB.ReadGeno.R --GenoFile /data1/UK_Biobank/ukb23155_b0_v1_s200604/ukb23155_c10_b0_v1.bed --OutputPrefix output --ControlFile ControlFile.txt

#### R ####
# IDsToIncludeFile = "/home/wenjianb/Docker/GRAB/IDsToIncludeFile.txt"
# ControlFile = "/home/wenjianb/Docker/GRAB/ControlFile.txt"
# write.table(c("10:133626803:I:1","10:133626804:T:C"), IDsToIncludeFile, row.names=F, col.names=F, quote=F, sep="\t")
# control.dataframe = data.frame(name = "IDsToIncludeFile",
#                                value = IDsToIncludeFile)
# write.table(control.dataframe, ControlFile, col.names=F, row.names=F, quote=F, sep="\t")

options(stringsAsFactors=F)

## load R libraries
library(GRAB)
library(optparse)

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option("--GenoFile", type="character", default="",
              help="a character of genotype file, bed file or bgen file"),
  make_option("--GenoMarkerIndexFile", type="character", default=NULL,
              help="marker index file corresponding to GenoFile, bim file or bgi file."),
  make_option("--GenoSampleFile", type="character", default=NULL,
              help="sample file corresponding to GenoFile, fam file or sample file."),
  make_option("--SampleIDFile", type="character", default=NULL,
              help="a file to store sample IDs, one column (no header)."),
  make_option("--IDsToIncludeFile", type="character", default=NULL,
              help="a file to store control list, two tab-seperated columns (no headers): column 1 is name and column 2 is value."),
  make_option("--OutputPrefix", type="character", default=NULL,
              help="a file to store control list, two tab-seperated columns (no headers): column 1 is name and column 2 is value.")
)


## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

extractSampleIDs = function(SampleIDFile)
{
  if(is.null(SampleIDFile))
    return(NULL);
  
  SampleIDs = read.table(SampleIDFile, header = F)
  if(ncol(SampleIDs) != 1)
    stop("ncol(SampleIDs) != 1")
  
  SampleIDs = as.character(unlist(SampleIDs))
  return(SampleIDs)
}

# extractControl = function(ControlFile)
# {
#   if(is.null(ControlFile))
#     return(NULL)
#   
#   control.dataframe = read.table(ControlFile, header = F, sep = "\t")
#   if(ncol(control.dataframe) != 2)
#     stop("ncol(control.dataframe) != 2")
#   
#   m = nrow(control.dataframe)
#   if(m == 0)
#     stop("nrow(control.dataframe) == 0")
#   
#   control.list = list()
#   for(i in 1:m)
#   {
#     name = control.dataframe[i, 1]
#     value = control.dataframe[i, 2]
#     control.list[[name]] = value
#   }
#   
#   return(control.list)
# }

SampleIDs = extractSampleIDs(opt$SampleIDFile)
# control = extractControl(opt$ControlFile)

if(is.null(opt$OutputPrefix))
  stop("argument of OutputPrefix is required.")

OutputFile.GenoMat = paste0(opt$OutputPrefix, ".GenoMat.txt")
OutputFile.markerInfo = paste0(opt$OutputPrefix, ".markerInfo.txt")

GenoFileIndex = c(opt$GenoMarkerIndexFile, opt$GenoSampleFile)

GenoList = GRAB.ReadGeno(opt$GenoFile,
                         GenoFileIndex,
                         SampleIDs,
                         control = list(IDsToIncludeFile = opt$IDsToIncludeFile),
                         sparse = FALSE)

GenoMat = data.table::data.table(GenoList$GenoMat, keep.rownames = T)
markerInfo = GenoList$markerInfo

# save(GenoList, 
#      file = paste0(opt$OutputPrefix,".RData"))

data.table::fwrite(GenoMat, OutputFile.GenoMat, 
                   col.names = T, row.names = F, sep = "\t", quote = F)
data.table::fwrite(markerInfo, OutputFile.markerInfo, 
                   col.names = T, row.names = F, sep = "\t", quote = F)
