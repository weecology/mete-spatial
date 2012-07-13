library(bigmemory)

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

S = round(10^seq(log10(10), log10(100), length.out=20))
N = round(10^seq(log10(120), log10(5e5), length.out=20))
B = 12

files = dir('./comms')

grains = 2^c(0, 2, 4, 6, 8)
ncomm = 200

for (s in S) {
  for (n in N) {
    fileSuffix = paste('_S', s, '_N', n, '_C200_B', B, '_grid', sep='')
    fileName = paste('simulated_comms', fileSuffix, '.txt', sep='')
    if (!(fileName %in% files)) 
      next
    comms = read.big.matrix(file.path('./comms', fileName), header=TRUE, 
                            type='integer', sep=',', descriptor = fileSuffix)
    Ns = n_pixels_long(B)
    Ms = n_pixels_wide(B)
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
    if (!exists('sarOut'))
      sarOut = data.frame(S = s, N=n, grains = grains, 
                        sr.lo = sarLo$richness, sr.avg = sarAvg$richness, 
                        sr.hi = sarHi$richness, ind.lo = sarLo$indiv, 
                        ind.avg = sarAvg$indiv, ind.hi = sarHi$indiv, 
                        count = sarAvg$count)
    else
      sarOut = rbind(sarOut,
                     data.frame(S = s, N = n, grains = grains, 
                                sr.lo = sarLo$richness, sr.avg = sarAvg$richness, 
                                sr.hi = sarHi$richness, ind.lo = sarLo$indiv, 
                                ind.avg = sarAvg$indiv, ind.hi = sarHi$indiv, 
                                count = sarAvg$count))
    write.csv(sarOut, file ='./sar/param_sar_avgs.csv', row.names=FALSE)
    print(paste('S =', s, 'N =', n, 'is done', sep=' '))
    rm(comms)
  }
}
write.csv(sarOut, file ='./sar/param_sar_avgs.csv', row.names=FALSE)

