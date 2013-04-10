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
  nperm = as.numeric(clArgs[4])
  npar = as.numeric(clArgs[5])
}
if (length(clArgs) == 0) {
  stop('Must specify commName, metricsToCalc, dataType, nperm, and npar at
       command line')
  q('no')
}
if (is.na(nperm)) {
  nperm = NULL
  npar = 1
}

## read in comms file 
comms = read.csv(paste('./data/',commName, '_comms.csv', sep=''))

## specify how to bin the spatial lags
spat_breaks = read.csv('./data/nbreaks.csv')
spat_breaks$nbreaks = spat_breaks$nbreaks + 1
if (any(spat_breaks$comm == commName)) {
  breaks = spat_breaks$nbreaks[spat_breaks$comm == commName]
  log = TRUE
}
if (!any(spat_breaks$comm == commName)) {
  breaks = NA
}

## specify quantiles to examine variogram at
quants = c(0.25, 0.50, 0.75)

## compute Dist Decay statistics
metrics = calcMetrics(comms, metricsToCalc, dataType, breaks=breaks, 
                      log=log, quants=quants, nperm=nperm, npar=npar,
                      writeToFile=TRUE, fileSuffix=commName)
