##Author: Dan McGlinn
##Purpose: to compare MaxEnt generate community to null models
##Outputs: the result of a null permutation routine.

library(vegan)
library(bigmemory)

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 1) {
  S = clArgs[1]
  N = clArgs[2]
  ncomm = clArgs[3]
  bisec_fine = as.numeric(clArgs[4])
  bisec_coarse = as.numeric(clArgs[5])
  grain_fine = as.numeric(clArgs[6]))
  transect = clArgs[7]
  dataType = clArgs[8]
  metricsToCalc = clArgs[9]
  direction = clArgs[10]
  tolerance = clArgs[11]
  name = clArgs[12]
  big = clArgs[13]
  nperm = clArgs[14]
  npar = clArgs[15]
}
if (!exists(as.character(substitute(S)))) {
  S = 10
  N = 100
  ncomm = 200
  bisec_fine = 12
  bisec_coarse = 6
  grain_fine = 1
  transect = FALSE
  dataType = 'abu'
  metricsToCalc = 'all'
  direction = 'NA'
  tolerance = 'NA'
  name = 'NA'
  big = FALSE
  nperm = NULL
  npar = 1
}

direction = ifelse(direction == 'NA', 'omnidirectional', as.numeric(direction))
tolerance = ifelse(tolerance == 'NA', NA, as.numeric(tolerance)) 

if(name == 'NA'){
  fileSuffix = ifelse(transect,
               paste('S', S, '_N', N, '_C', ncomm, '_B', bisec_fine, '_transect',
                     sep=''),
               paste('S', S, '_N', N, '_C', ncomm, '_B', bisec_fine, '_grid', sep=''))
}
if(name != 'NA'){
  fileSuffix = ifelse(transect,
               paste(name, '_C', ncomm, '_B', bisec_fine, '_transect', sep=''),
               paste(name, '_C', ncomm, '_B', bisec_fine, '_grid', sep=''))
}

fileName = paste('simulated_comms_', fileSuffix, '.txt', sep='')

big = ifelse(big, TRUE, FALSE)
if (big)
  comms = read.big.matrix(file.path('./comms', fileName), header=TRUE, 
                          type='integer', sep=',', descriptor = fileSuffix)
if (!big)
  comms = read.csv(file.path('./comms', fileName), header=T)

gc()

## specify how to bin the spatial lags
spat_breaks = read.csv('./data/nbreaks.csv')
spat_breaks$nbreaks = spat_breaks$nbreaks + 1
if (any(spat_breaks$comm == name)) {
  nbreaks = spat_breaks$nbreaks[spat_breaks$comm == name]
  breaks = sapply(nbreaks, function(x) list(c('log2', x)))
}
if (!any(spat_breaks$comm == name)) {
  breaks = NA
}

## specify quantiles to examine variogram at
quants = c(0.25, 0.50, 0.75)

## specify bisections levels
bisec = seq(bisec_fine, bisec_coarse, -2)

## specify grain names
grain_names = round(grain_fine * 2^(max(bisec) - bisec), 2)

## loop through all the communities
comm_ids = unique(comms[ , 1])
metrics = vector('list', length(comm_ids))
names(metrics) = comm_ids      
for (i in seq_along(comm_ids)) {
  comm_tmp = comms[comms[ , 1] == comm_ids[i], ]
  mat = comm_tmp[ , -(1:3)]
  coords = comm_tmp[ , 2:3]
  ## aggregate community matrix to other appropriate spatial lags
  comms_aggr = aggr_comm_matrix(mat, coords, bisec, grain_names)
  metrics[[i]] = calcMetrics(comms_aggr, metricsToCalc, dataType, breaks=breaks,
                        quants=quants, direction=direction, tolerance=tolerance,
                        nperm=nperm, npar=npar, writeToFile=FALSE)
  ## update save file as loop progresses
  save(metrics, file=paste('./', metricsToCalc, '/', metricsToCalc, '_',
                           fileSuffix, '_', dataType, '.Rdata', sep=''))
}


