## Author: Dan McGlinn
## Purpose: to analyze the spatial patterns of community structure
## in the empirical datasets

setwd('~/maxent/spat')

library(vegan)
source('./scripts/spat_sim_vario_func.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 1) {
  commName = clArgs[1]
  metricsToCalc = clArgs[2]
  dataType = clArgs[3]
  method = clArgs[4]
  npar = as.numeric(clArgs[5])
  nperm = as.numeric(clArgs[6])
}
if (length(clArgs) == 0) {
  stop('Must specify commName, metricsToCalc, dataType, nperm, and npar at
       command line')
  q('no')
}
if (method == 'uni') {
  univariate = TRUE
}
if (method == 'multi') {
  univariate = FALSE
}
if (npar > 1) {
  library(snowfall)
  library(rlecuyer)
  sfInit(parallel=TRUE, cpus=npar, type="SOCK")
}

## read in comms file 
comms = read.csv(paste('./data/',commName, '_comms.csv', sep=''))

## specify quantiles to examine variogram at
quants = c(0.25, 0.50, 0.75)

## compute Dist Decay statistics
metrics = calc_metrics_bisect(comms, metricsToCalc, dataType, quants,
                              nperm, univariate, writeToFile=TRUE,
                              fileSuffix=commName)

sfStop()
