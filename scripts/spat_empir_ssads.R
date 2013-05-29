setwd('~/maxent/spat')
source('./scripts/spat_functions.R')

## read in empirical data---------------------------------
fileNames = dir('./data')
commFiles = grep('comms', fileNames)
dat = vector('list', length(commFiles))
names(dat) = sub('_comms.csv', '', fileNames[commFiles])
for (i in seq_along(commFiles))
  dat[[i]] = as.matrix(read.csv(paste('./data/', fileNames[commFiles[i]], sep='')))

## compute SSAD at each spatial scale---------------------
ssads = sapply(dat, get_SSAD)

## export SSADs as .csv files-----------------------------
for (i in seq_along(dat))
  write.csv(ssads[[i]], 
            file=paste('./data/', names(dat)[i], '_ssads.csv',sep=''), 
            row.names=F)


  