
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
  make_option("--objNullFile", type="character", default="",
              help="a file of object to load the results in GRAB.NullModel."),
  make_option("--GenoFile", type="character", default="",
              help="a character of genotype file, bed file or bgen file"),
  make_option("--GenoMarkerIndexFile", type="character", default=NULL,
              help="marker index file corresponding to GenoFile, bim file or bgi file."),
  make_option("--GenoSampleFile", type="character", default=NULL,
              help="sample file corresponding to GenoFile, fam file or sample file."),
  make_option("--GroupFile", type="character", default=NULL,
              help="a file to store group information for all regions."),
  make_option("--SparseGRMFile", type="character", default=NULL,
              help="a file to store sparse GRM information, three tab-seperated columns (no headers): columns 1 and 2 are subject IDs and column 3 is the kinship coefficient."),
  make_option("--ControlFile", type="character", default=NULL,
              help="a file to store information in control list, two tab-seperated columns (no headers): columns 1 is the list key and column 2 is the list value."),
  make_option("--MaxMAFVec", type="character", default="0.05,0.01,0.001,0.0001",
              help="max MAF cutoff to include markers for region-level analysis. Multiple cutoffs separated with ';' is supported."),
  make_option("--annoVec", type="character", default="lof,lof:missense,lof:missense:synonymous",
              help="annotation list used in analysis. Multiple annotations seperated with ';' is supported"),
  make_option("--OutputPrefix", type="character", default=NULL,
              help="a file to store control list, two tab-seperated columns (no headers): column 1 is name and column 2 is value.")
)


## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

### output files

# OutputFile
# paste0(OutputFile, ".index");
# paste0(OutputFile, ".markerInfo")
# paste0(OutputFile, ".otherMarkerInfo")

###############

load(opt$objNullFile)

if(!"objNull" %in% ls())
  stop('!"objNull" %in% ls()', ":\t objNullFile include includes an R object named with 'objNull'.")

if(is.null(opt$OutputPrefix))
  stop("argument of OutputPrefix is required.")

GenoFileIndex = c(opt$GenoMarkerIndexFile, opt$GenoSampleFile)
OutputFile = opt$OutputPrefix;

extractControl = function(ControlFile)
{
  if(is.null(ControlFile))
    return(NULL)
  
  control.dataframe = read.table(ControlFile, header = F, sep = "\t")
  if(ncol(control.dataframe) != 2)
    stop("ncol(control.dataframe) != 2")
  
  m = nrow(control.dataframe)
  if(m == 0)
    stop("nrow(control.dataframe) == 0")
  
  control.list = list()
  for(i in 1:m)
  {
    name = control.dataframe[i, 1]
    value = control.dataframe[i, 2]
    control.list[[name]] = value
  }
  
  return(control.list)
}

###############

# MaxMAFVec = as.numeric(strsplit(opt$MaxMAF, split=",")[[1]])
# if(any(is.na(MaxMAFVec)))
#   stop("any(is.na(MaxMAFVec)): please check argument --MaxMAF carefully.")
# 
# annoVec = strsplit(opt$anno, split = ",")[[1]]
control = extractControl(opt$ControlFile)

GRAB.Region(objNull,
            opt$GenoFile,
            GenoFileIndex,
            OutputFile,
            OutputFileIndex = NULL,
            opt$GroupFile,
            opt$SparseGRMFile,
            MaxMAFVec = opt$MaxMAFVec,  
            annoVec = opt$annoVec,    
            control = control)
