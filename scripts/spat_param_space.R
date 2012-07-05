options(scipen=10) 

setwd('~/maxent/spat/')
source('./scripts/spat_sim_vario_func.R')

S = round(10^seq(log10(10), log10(100), length.out=20))
N = round(10^seq(log10(120), log10(5e5), length.out=20))

simSorAbu = loadSimResults(S, N, './sorensen/')
simSorAbuAvg = avgSimResults(simSorAbu, 'sorensen')

#save(simSorAbuAvg, file='./sorensen/simSorAbuAvg.Rdata')
#load('./sorensen/simSorAbuAvg.Rdata')

stats = getSimStats(simSorAbuAvg, S, N)

pdf('./figs/parm_space_r2_sorensen_abu.pdf', width = 7 * 2, height = 7)
  lims = range(as.vector(stats[ ,'r2', , , , ]),na.rm=TRUE)
  nbrks = 12
  brks = seq(lims[1], lims[2], length.out=nbrks)
  par(mfrow=c(1, 2))
  for (meth in c('ols', 'wtr')) { 
    hist(stats['exp', 'r2', meth, , , ], breaks = brks, ylim=c(0, 1.7e3), main=meth,
         col='red', xlab=expression(R^2 * ' of model'))
    hist(stats['pwr', 'r2', meth, , , ], breaks = brks, ylim=c(0, 1.7e3), 
         add=TRUE, col='blue')
    legend('topleft', c('Exponential Model', 'Power Model'), col=c('red', 'blue'),
           lwd=8, cex=2, bty='n')
  }
dev.off()


pwrStats = drop(stats['pwr', , , , ,])

## plot intercept,slope,and R2 for pwr model of DD at each grain

pdf('./figs/stats_sor_abu_all_grains_pwr_wtr_.pdf', width=7 * 2, height=7)
  for (g in seq_along(grains)) {
    par(mfrow=c(1, 3))
    for (i in 1:3) {
      image(pwrStats[i, 'wtr', g, , ])
    }  
  }
dev.off()

pdf('./figs/stats_sor_abu_grain1_pwr_wtr_.pdf', width=7 * 2, height=7)
  par(mfrow=c(1, 3))
  for (i in 1:3) {
    image(pwrStats[i, 'wtr', 1, , ])
  }  
dev.off()

##create legend figure
par(mfrow=c(1,3))
for(i in 1:3){
  tmp = as.vector(pwrStats[i, 'ols', 1, , ])
  x = round(seq(min(tmp,na.rm=T), max(tmp,na.rm=T), length.out=7), 3)
  image(matrix(x),axes=F)
  axis(side=1,at=seq(0,1,length.out=7),labels=x,tick=F,cex.axis=2)
}

svals = rep(S, length(N))
nvals = rep(N, each=length(S))

pdf('./figs/log_n_over_s_coef_sor_abu_grain1_pwr_wtr.pdf', width = 7 * 2, height = 7)
  ratio = log(nvals / svals)
  ratio = array(ratio, dim=c(length(S), length(N)))
  xlab = 'log(N/S)'
  g = 1
  meth = 'wtr'
  lwd = 3
  par(mfrow=c(2,3))
  plot(ratio, pwrStats['b0', meth, g, , ], xlab=xlab, ylab='Intercept', pch=19)
  plot(ratio, pwrStats['b0', meth, g, , ], xlab=xlab, ylab='Intercept', type='n')
  for (s in seq_along(S))
    lines(ratio[s, ], pwrStats['b0', meth, g, s, ], col='palevioletred', lwd=lwd)
  plot(ratio, pwrStats['b0', meth, g, , ], xlab=xlab, ylab='Intercept', type='n')
  for (n in seq_along(N))
    lines(ratio[ , n], pwrStats['b0', meth, g, , n], col='dodgerblue', lwd=lwd)
  ##
  plot(ratio, pwrStats['b1', meth, g, , ], xlab=xlab, ylab='Slope', pch=19)
  plot(ratio, pwrStats['b1', meth, g, , ], xlab=xlab, ylab='Slope', type='n')
  for (s in seq_along(S))
    lines(ratio[s, ], pwrStats['b1', meth, g, s, ], col='palevioletred', lwd=lwd)
  plot(ratio, pwrStats['b1', meth, g, , ], xlab=xlab, ylab='Slope', type='n')
  for (n in seq_along(N))
    lines(ratio[ , n], pwrStats['b1', meth, g, , n], col='dodgerblue', lwd=lwd)
dev.off()

## consider 3 of the grains
pdf('./figs/log_n_over_s_coef_sor_abu_3grains_pwr_wtr.pdf', width = 7 * 2,
    height = 7 * 2)
  ratio = log(nvals / svals)
  ratio = array(ratio, dim=c(length(S), length(N)))
  xlab = 'log(N/S)'
  par(mfrow=c(3, 3))
  for (cof in c('b0', 'b1')) {
    for (g in c(1, 3, 5)) {
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, pch=19,
           main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(s in seq_along(S))
        lines(ratio[s, ], pwrStats[cof, meth, g, s, ], col='palevioletred', 
              lwd=lwd)
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(n in seq_along(N))
        lines(ratio[ , n],pwrStats[cof, meth, g, , n], col='dodgerblue', lwd=lwd)
    }  
  }  
dev.off()

## overlap the grains onto a single panel
pdf('./figs/log_n_over_s_coef_sor_abu_all_grains_overlap_pwr_wtr.pdf', 
    width = 7 * 2, height = 7)
  ratio = log(nvals / svals)
  ratio = array(ratio, dim=c(length(S), length(N)))
  xlab = 'log(N/S)'
  par(mfrow=c(1, 2))
  for (cof in c('b0', 'b1')) {
    plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n',
         ylim = range(pwrStats[cof, meth, , , ], na.rm=T))
      for (g in seq_along(grains)) {
        points(ratio, pwrStats[cof, meth, g, , ], col=g, pch=1)
      }  
    if (cof == 'b0')
      legend('bottomright', legend=c('Grains', grains), col=c(NA, 1:length(grains)),
             pch=1, bty='n', cex = 2, lwd=3, lty=NA)
  }  
dev.off()



## consider log log
pdf('./figs/log_log_n_over_s_coef_sor_abu_3grains_pwr_wtr.pdf', width = 7 * 2,
    height = 7 * 2)
  ratio = log(log(nvals / svals))
  ratio = array(ratio, dim=c(length(S), length(N)))
  xlab = 'log(log(N/S))'
  par(mfrow=c(3, 3))
  for (cof in c('b0', 'b1')) {
    for (g in c(1, 3, 5)) {
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, pch=19,
           main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(s in seq_along(S))
        lines(ratio[s, ], pwrStats[cof, meth, g, s, ], col='palevioletred', 
              lwd=lwd)
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(n in seq_along(N))
        lines(ratio[ , n],pwrStats[cof, meth, g, , n], col='dodgerblue', lwd=lwd)
    }  
  }  
dev.off()

## log transforms of N

ratio = log(nvals)
ratio = array(ratio, dim=c(length(S), length(N)))
xlab = 'log(N)'
pdf('./figs/logN_coef_sor_abu_3grains_pwr_wtr.pdf', width = 7 * 2, height = 7 * 2)
  for (cof in c('b0', 'b1')) {
    par(mfrow=c(3, 3))
    for (g in c(1, 3, 5)) {
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, pch=19,
           main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(s in seq_along(S))
        lines(ratio[s, ], pwrStats[cof, meth, g, s, ], col='palevioletred', 
              lwd=lwd)
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(n in seq_along(N))
        lines(ratio[ , n],pwrStats[cof, meth, g, , n], col='dodgerblue', lwd=lwd)
    }  
  }  
dev.off()

## log transforms of S

ratio = log(svals)
ratio = array(ratio, dim=c(length(S), length(N)))
xlab = 'log(S)'
pdf('./figs/logS_coef_sor_abu_3grains_pwr_wtr.pdf', width = 7 * 2, height = 7 * 2)
  for (cof in c('b0', 'b1')) {
    par(mfrow=c(3, 3))
    for (g in c(1, 3, 5)) {
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, pch=19,
           main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(s in seq_along(S))
        lines(ratio[s, ], pwrStats[cof, meth, g, s, ], col='palevioletred', 
              lwd=lwd)
      plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
      for(n in seq_along(N))
        lines(ratio[ , n],pwrStats[cof, meth, g, , n], col='dodgerblue', lwd=lwd)
    }  
  }  
dev.off()

####
## create flat file to do visual checks that other plots
## are correct
meth = 'wtr'
g = 1
flat = data.frame(S = svals, N = nvals, ratio = as.vector(ratio), 
                  b0 = as.vector(pwrStats['b0', meth, g, , ]),
                  b1 = as.vector(pwrStats['b1', meth, g, , ]))

par(mfrow=c(1,2))
plot(b0 ~ ratio, data=flat, type='n')
lines(b0 ~ ratio, data=flat, subset= S == 10)
plot(b0 ~ ratio, data=flat, type='n')
lines(b0 ~ ratio, data=flat, subset= N == 120)
## everytning matches with what the plotting routines below produce


