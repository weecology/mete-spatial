## Purpose: to compute the empirical SARs

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

print('Computing empirical SARs...')

fileNames = dir('./data')
commFiles = grep('comms', fileNames)

## First aggregate community files into a single list
dat = vector('list', length(commFiles))
names(dat) = sub('_comms.csv', '', fileNames[commFiles])
for (i in seq_along(commFiles))
  dat[[i]] = as.matrix(read.csv(paste('./data/', fileNames[commFiles[i]], sep='')))

Amin = unlist(lapply(dat, function(x) unique(x[ , 1])[1]))
for (i in seq_along(dat))
  dat[[i]] = (dat[[i]][dat[[i]][ , 1] == Amin[i], ])

## load empirical constants
shrtnames = read.table('./data/shrtnames.txt', colClasses='character')
grain_fine = read.table('./data/grain_fine.txt')
bisect_fine = read.table('./data/bisect_fine.txt')
bisect_coarse = read.table('./data/bisect_coarse.txt')

proper_order = match(names(dat), shrtnames)
grain_fine = grain_fine[proper_order]
bisect_fine = bisect_fine[proper_order]
bisect_coarse = bisect_coarse[proper_order]

## convert them to multidimensional arrays

Ns = n_pixels_long(bisect_fine)
Ms = n_pixels_wide(bisect_fine)
psp = vector('list',length(dat))
names(psp) = names(dat)
for(i in seq_along(dat))
  psp[[i]] = mat2psp(dat[[i]][ , -(1:3)], dat[[i]][ , 2:3],
                     Ns[i], Ms[i])

## compute SARs, non-movinging window 
grains = sapply(1:length(psp), function(x) 2^(0:bisect_fine[1, x]))

sar = vector('list',length(psp))
names(sar) = names(psp)
for (i in seq_along(psp)) {
  sar[[i]] = getSAR(psp[[i]], grains[[i]])
  ## add area m2 column
  sar[[i]] = data.frame(sar[[i]], area = sar[[i]][ , 1] * grain_fine[1, i])
}

## export results
for (i in seq_along(sar)) 
  write.csv(sar[[i]], file=paste('./sar/', names(sar)[i], '_empir_sar.csv', 
                                sep=''), row.names=FALSE)

## create flat file of sar data
for (i in seq_along(sar)) {
  if (i == 1)
    dat = data.frame(site=names(sar)[i], sar[[i]])
  else
    dat = rbind(dat, data.frame(site=names(sar)[i], sar[[i]]))
}

write.csv(dat, file='./sar/empir_sars.csv', row.names=FALSE)

print('Computing empirical SARs, complete!')
