## Author: Dan McGlinn
## Purpose: to analyze the spatial patterns of community structure
## in the sherman census 3 which was surveyed in 1999 for each spatial
## scale of interest
## Metadata: sherman_files.doc
setwd('/home/danmcglinn/maxent/spat')
.libPaths('/home/danmcglinn/R/x86_64-pc-linux-gnu-library/2.14')

library(vegan)
library(danspkg)
source('spat_sim_vario_func.R')

clArgs <- commandArgs(TRUE)
if( length(clArgs) > 1){
  metricsToCalc <- clArgs[1]
  dataType <- clArgs[2]
}
if( !exists(as.character(substitute(metricsToCalc))) ){
  metricsToCalc <- 'all'
  dataType <- 'abu'
}

## read in comms file for the 1999 census (i.e., census3) of sherman
comms = read.csv('./data/sherman_comms.csv',header=TRUE)

## compute  Dist Decay statistics
metrics = calcMetrics(comms,metricsToCalc,dataType,writeToFile=TRUE,
                      fileSuffix='_sherman')






