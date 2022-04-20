#!/usr/bin/Rscript

options(stringsAsFactors=F)

install.packages("optparse", repos = "https://cloud.r-project.org")

library(optparse)

## set list of cmd line arguments
option_list <- list(
  make_option("--GenoFile", type = "character", default = "",
    help="a character of genotype file."),
  make_option("--GenoFileIndex", type="character", default="",
    help="additional index file(s) corresponding to GenoFile."),
  make_option("--control-IDsToIncludeFile", type="character", default="",
    help="a file of marker IDs to include, one column (no header)."))

opt <- parse_args(OptionParser(option_list=option_list))

# print some progress messages to stderr if "quietly" wasn't requested

print(opt)

# ./grab.R --GenoFile=123 --control-IDsToIncludeFile xxx --control-IDsToIncludeFILE yyy
# Error in getopt_options(object, args) :
#   Error in getopt(spec = spec, opt = args) :
#   long flag "control-IDsToIncludeFILE" is invalid
# Calls: parse_args -> parse_options -> getopt_options


# [wenjianb@master Docker]$ ./grab.R --GenoFile=123 --control-IDsToIncludeFile xxx
# $GenoFile
# [1] "123"
# 
# $GenoFileIndex
# [1] ""
# 
# $`control-IDsToIncludeFile`
# [1] "xxx"
# 
# $help
# [1] FALSE

opt = list(GenoFile = "xx",
           `control-arg1` = "xx",
           `control-arg2` = "xx")

opt.names = names(opt)

opt.control = grep("control-", opt.names, value = T)
# opt.control = gsub("control-", "", opt.control)

control.list = opt[opt.control]
control.list
names(control.list) = gsub("control-", "", names(control.list))
control.list

