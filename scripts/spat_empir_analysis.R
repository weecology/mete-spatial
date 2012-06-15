## Author: Dan McGlinn
## Purpose: to analyze the spatial patterns of community structure
## in the empirical datasets

setwd('~/maxent/spat')

library(vegan)
library(danspkg)
source('./scripts/spat_sim_vario_func.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 1) {
  commName = clArgs[1]
  metricsToCalc = clArgs[2]
  dataType = clArgs[3]
}
if (length(clArgs) == 0) {
  stop('Must specify commName, metricsToCalc, & dataType at command line')
  q('no')
}

## read in comms file 
comms = read.csv(paste('./data/',commName, '_comms.csv', sep=''))

## specify quantiles to examine variogram at
quants = c(0.25, 0.50, 0.75)

## compute Dist Decay statistics
metrics = calcMetrics(comms, metricsToCalc, dataType, quants=quants,
                      writeToFile=TRUE, fileSuffix=commName)






