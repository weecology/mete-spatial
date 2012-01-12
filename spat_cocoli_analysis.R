## Author: Dan McGlinn
## Purpose: to analyze the spatial patterns of community structure
## in the cocoli census 3 which was surveyed in 1998 for each spatial
## scale of interest
## Metadata: cocoli_files.doc
setwd('/home/danmcglinn/maxent/spat')

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

## read in comms file for the 1998 census (i.e., census3) of cocoli
comms = read.csv('./data/cocoli_comms.csv',header=TRUE)

## compute  Dist Decay statistics
metrics = calcMetrics(comms,metricsToCalc,dataType,writeToFile=TRUE,
                      fileSuffix='_cocoli')






