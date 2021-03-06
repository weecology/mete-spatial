##Purpose: to compute the DDR of the METE simulated community matrices

print('METE simulated DDR analysis, ...')

library(vegan)
library(bigmemory)

source('./scripts/spat_functions.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 1) {
  S = clArgs[1]
  N = clArgs[2]
  ncomm = clArgs[3]
  bisect_fine = as.numeric(clArgs[4])
  bisect_coarse = as.numeric(clArgs[5])
  grain_fine = as.numeric(clArgs[6])
  transect = clArgs[7]
  dataType = clArgs[8]
  metricsToCalc = clArgs[9]
  bisect = clArgs[10]
  direction = clArgs[11]
  tolerance = clArgs[12]
  name = clArgs[13]
  big = ifelse(clArgs[14], TRUE, FALSE)
  iteration = as.numeric(clArgs[15])
}
if (!exists(as.character(substitute(S)))) {
  S = 10
  N = 100
  ncomm = 200
  bisect_fine = 12
  bisect_coarse = 6
  grain_fine = 1
  transect = FALSE
  dataType = 'abu'
  metricsToCalc = 'all'
  bisect = FALSE
  direction = 'NA'
  tolerance = 'NA'
  name = 'NA'
  big = FALSE
  iteration = 1
}

direction = ifelse(direction == 'NA', 'omnidirectional', as.numeric(direction))
tolerance = ifelse(tolerance == 'NA', NA, as.numeric(tolerance)) 

if(name == 'NA'){
  fileSuffix = ifelse(transect,
               paste('S', S, '_N', N, '_C', ncomm, '_B', bisect_fine, '_transect',
                     sep=''),
               paste('S', S, '_N', N, '_C', ncomm, '_B', bisect_fine, '_grid', sep=''))
}
if(name != 'NA'){
  fileSuffix = ifelse(transect,
               paste(name, '_C', ncomm, '_B', bisect_fine, '_transect', sep=''),
               paste(name, '_C', ncomm, '_B', bisect_fine, '_grid', sep=''))
}


fileName = paste('simulated_comms_', fileSuffix, '.txt', sep='')

if (big)
  comms = read.big.matrix(file.path('./comms', fileName), header=TRUE, 
                          type='short', sep=',', descriptor = fileSuffix)
if (!big)
  comms = read.csv(file.path('./comms', fileName), header=T)

gc()

if(bisect) {
  fileSuffix = paste(fileSuffix, '_bisect', sep='')
}

## specify quantiles to examine variogram at
quants = c(0.25, 0.50, 0.75)

## specify bisections levels
bisect_lvs = seq(bisect_fine, bisect_coarse, -2)

## specify grain names
grain_names = round(grain_fine * 2^(max(bisect_lvs) - bisect_lvs), 2)

## loop through all the communities
comm_ids = unique(comms[ , 1])

if (iteration > 1) {
  load(paste('./', metricsToCalc, '/', metricsToCalc, '_',
             fileSuffix, '_', dataType, '.Rdata', sep=''))
}
if (iteration == 1) {
  metrics = vector('list', length(comm_ids))
  names(metrics) = comm_ids      
}

for (i in seq_along(comm_ids)) {
  if (i >= iteration) {
    comm_tmp = comms[comms[ , 1] == comm_ids[i], ]
    mat = comm_tmp[ , -(1:3)]
    coords = comm_tmp[ , 2:3]
    ## aggregate community matrix to other appropriate spatial lags
    comms_aggr = aggr_comm_matrix(mat, coords, bisect_lvs, grain_names)
    if (bisect)
      metrics[[i]] = calc_metrics_bisect(comms_aggr, metricsToCalc, dataType, 
                                         quants, writeToFile=FALSE)
    else
      metrics[[i]] = calc_metrics(comms_aggr, metricsToCalc, dataType, breaks=breaks,
                                  log=log, quants=quants, direction=direction,
                                  tolerance=tolerance,  writeToFile=FALSE)
    ## update save file as loop progresses
    save(metrics, file=paste('./', metricsToCalc, '/', metricsToCalc, '_',
                             fileSuffix, '_', dataType, '.Rdata', sep=''))
  }
}

print('METE simulated DDR analysis, complete!')

