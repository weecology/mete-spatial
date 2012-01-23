## Author: Dan McGlinn
## Purpose: to analyze the spatial patterns of community structure
## in the BCI census 7 which was surveyed in 2010 for each spatial
## scale of interest
## Metadata: bci50ha.doc
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

## read in comms file for the 2010 census (i.e., census7) of BCI
comms = read.csv('./data/bci_comms.csv',header=TRUE)

## compute  Dist Decay statistics
metrics = calcMetrics(comms,metricsToCalc,dataType,writeToFile=TRUE,fileSuffix='_bci')






