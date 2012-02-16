## Purpose: to compute the empirical SARs

setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

fileNames = dir('./sar')
meteFiles = grep('mete',fileNames)
mete = vector('list',length(meteFiles))
for(i in seq_along(meteFiles))
  mete[[i]] = read.csv(paste('./sar/',fileNames[meteFiles[i]],sep=''))
names(mete) = sub('_mete_sar.txt','',fileNames[meteFiles])

##for cocoli1, cocoli2, sherman1, sherman2 we need to adjust the sars
##b/c the mete prediction fails at scales where there is less than 1 individual
toFix = c('cocoli1', 'cocoli2', 'sherman1', 'sherman2')
toFix = which(names(mete)%in%toFix)
for(i in seq_along(toFix)){
  mete[[toFix[i]]]$area = mete[[toFix[i]]]$area * 2
  mete[[toFix[i]]] = rbind(c(1,NA),mete[[toFix[i]]])
}
toFix = 'cross'
toFix = which(names(mete)%in%toFix)
for(i in seq_along(toFix)){
  mete[[toFix[i]]]$area = mete[[toFix[i]]]$area * 4
  mete[[toFix[i]]] = rbind(c(1,NA),c(2,NA),mete[[toFix[i]]])
}

fileNames = dir('./data')
commFiles = grep('comms',fileNames)
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
Ns = c(rep(128,3),64,16,rep(128,2),64)
Ms = c(rep(64,4),16,rep(64,2),64)
psp = vector('list',length(dat))
for(i in seq_along(dat))
  psp[[i]] = mat2psp(dat[[i]][,-(1:3)],Ns[i],Ms[i])
names(psp) = names(dat)

AminExact = c(1e3/128 * 5e2/64,
              2e2/128 * 1e2/64,
              2e2/128 * 1e2/64,
              2e2/64 * 2e2/64,
              1,
              2e2/128 * 1e2/64,
              2e2/128 * 1e2/64,
              (1.4e2/64)^2)
## compute SARs, non-movinging window 
sar = vector('list',length(psp))
for(i in seq_along(psp)){
  sar[[i]] = getSAR(psp[[i]],mete[[i]]$area)
  ## scale area appropriately
  sar[[i]][,1] = sar[[i]][,1]*AminExact[i]
}
names(sar) = names(psp)

## export results
for(i in seq_along(sar)){
  write.csv(sar[[i]],file=paste('./sar/',names(sar)[i],'_empir_sar.csv',sep=''),
            row.names=FALSE)
}



