## Purpose: to compute the empirical SARs

setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

fileNames = dir('./data')
commFiles = grep('comms',fileNames)
commFiles = commFiles[-c(4,9)] ## needed for the time being until we clear out the old comm files
## First aggregate community files into a single list
dat = vector('list',length(commFiles))
for(i in seq_along(commFiles))
  dat[[i]] = read.csv(paste('./data/',fileNames[commFiles[i]],sep=''))
names(dat) = sub('_comms.csv','',fileNames[commFiles])

## drop abundance and data from the community files for grains larger than the minimum grain
Amin = unlist(lapply(dat,function(x)unique(x$grain)[1]))
for(i in seq_along(dat))
  dat[[i]] = (dat[[i]][dat[[i]]$grain == Amin[i],] > 0) * 1
gc()

## convert them to multidimensional arrays
unlist(lapply(dat,nrow))
Ns = c(rep(128,3),16,rep(128,2),64)
Ms = c(rep(64,3),16,rep(64,2),64)
psp = vector('list',length(dat))
for(i in seq_along(dat))
  psp[[i]] = mat2psp(dat[[i]][,-(1:3)],Ns[i],Ms[i])
names(psp) = names(dat)

## compute SARs
sar = vector('list',length(psp))
for(i in seq_along(psp))
  sar[[i]] = getSAR(psp[[i]],mete[[i]]$area)
names(sar) = names(psp)

## export results
for(i in seq_along(sar)){
  write.csv(sar[[i]],file=paste('./sar/',names(sar)[i],'_empir_sar.csv',sep=''),
            row.names=FALSE)
}

