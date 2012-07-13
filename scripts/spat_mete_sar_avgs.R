library(bigmemory)

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

fileNames = dir('./sar')

meteFiles = grep('mete_sar.txt', fileNames)
mete = vector('list', length(meteFiles))
names(mete) = sub('_mete_sar.txt', '', fileNames[meteFiles])
for( i in seq_along(meteFiles)) {
  mete[[i]] = read.csv(paste('./sar/', fileNames[meteFiles[i]], sep=''))
}

bisect_fine = read.table('./data/bisect_fine.txt')
shrtnames = read.table('./data/shrtnames.txt', colClasses='character')

files = dir('./comms')
ncomm = 200

for (i in seq_along(mete)) {
  name = names(mete)[i]
  bisect = as.numeric(bisect_fine[match(name, shrtnames)])
  fileSuffix = paste(name, '_C', ncomm, '_B', bisect, '_grid', sep='')
  fileName = paste('simulated_comms_', fileSuffix, '.txt', sep='')
  if (!(fileName %in% files)) 
    next
  comms = read.big.matrix(file.path('./comms', fileName), header=TRUE, 
                          type='integer', sep=',', descriptor = fileSuffix)
  Ns = n_pixels_long(bisect)
  Ms = n_pixels_wide(bisect)
  grains = mete[[i]]$area
  for (j in 1:ncomm) {
    comm_tmp = comms[comms[ , 1] == j, ]
    psp = mat2psp(comm_tmp[ , -(1:3)], Ns, Ms)
    if (j == 1) 
      sar = data.frame(comm = j, getSAR(psp, grains))
    else
      sar = rbind(sar, data.frame(comm = j, getSAR(psp, grains)))
  }  
  sarAvg = aggregate(cbind(richness, indiv, count) ~ grains, data = sar,
                     FUN = mean)
  sarHi = aggregate(cbind(richness, indiv, count) ~ grains, data = sar,
                    FUN = quantile, 0.975)
  sarLo = aggregate(cbind(richness, indiv, count) ~ grains, data = sar,
                    FUN = quantile, 0.025)
  out = data.frame(comm = name, grains = grains, 
                   sr.lo = sarLo$richness, sr.avg = sarAvg$richness, 
                   sr.hi = sarHi$richness, ind.lo = sarLo$indiv, 
                   ind.avg = sarAvg$indiv, ind.hi = sarHi$indiv, 
                   count = sarAvg$count)
  write.csv(out, file=paste('./sar/', name, '_mete_sar_avgs.csv', sep=''),
            row.names=FALSE)
  if (i == 1)
    sarOut = out
  else
    sarOut = rbind(sarOut, out)
  print(paste(name, 'is done', sep=' '))
  rm(comms)
  gc()
}

write.csv(sarOut, file ='./sar/mete_sar_avgs.csv', row.names=FALSE)

