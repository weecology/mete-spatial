options(scipen=10) 

setwd('~/maxent/spat/')
source('./scripts/spat_sim_vario_func.R')

S = round(10^seq(log10(10), log10(100), length.out=20))
N = round(10^seq(log10(120), log10(5e5), length.out=20))

#simSorAbu = loadSimResults(S, N, './sorensen/')
#simSorAbuAvg = avgSimResults(simSorAbu, 'sorensen')
#save(simSorAbuAvg, file='./sorensen/simSorAbuAvg.Rdata')

load('./sorensen/simSorAbuAvg.Rdata')

stats = getSimStats(simSorAbuAvg, S, N)
pwrStats = drop(stats['pwr', , , , ,])

grains = unique(simSorAbuAvg[[1]]$grain)
svals = rep(S, length(N))
nvals = rep(N, each=length(S))

ddr = read.csv('./sar/param_sar_avgs.csv')
## add the stats information to this flat file
meth = 'wtr'
for (g in seq_along(grains)) {
  for (s in seq_along(S)) {
    for (n in seq_along(N)) {
      true = ddr$grains == grains[g] & ddr$S == S[s] & ddr$N == N[n]
      ddr$b0[true] = pwrStats['b0', meth, g, s, n]
      ddr$b1[true] = pwrStats['b1', meth, g, s, n]
      ddr$r2[true] = pwrStats['r2', meth, g, s, n]
    }
  }
} 

## write out flat file
write.csv(ddr, file='./sorensen/param_ddr_wtr_pwr_stats.csv',
          row.names = FALSE)

##----------------------------------------------------------------------------



## create an example graphic for the simulated pattern
S[20]
N[20]
index = 20 * 20
tmp = simSorAbuAvg[[index]]
plot(avg ~ Dist, data=tmp, log='xy', type='n', frame.plot=F, axes=F, xlab='', ylab='')
axis(side=1, cex.axis=1.5, lwd=2)
axis(side=2, cex.axis=1.5, lwd=2)
for (g in seq_along(grains)) {
  points(avg ~ Dist, data=tmp, subset=grains==grains[g])
  abline(a=stats['pwr', 'b0', 'wtr', g, 20, 20], 
         b=stats['pwr', 'b1', 'wtr', g, 20, 20])
}

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
  meth = 'wtr'
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
  g = 1
  meth = 'wtr'
  xlab = 'log(N/S)'
  par(mfrow=c(1, 2))
  for (cof in c('b0', 'b1')) {
    plot(ratio , pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n',
         ylim = range(pwrStats[cof, meth, , , ], na.rm=T))
      for (g in seq_along(grains)) {
        points(ratio , pwrStats[cof, meth, g, , ], col=g, pch=1)
      }  
    if (cof == 'b0')
      legend('bottomright', legend=c('Grains', grains), col=c(NA, 1:length(grains)),
             pch=1, bty='n', cex = 2, lwd=3, lty=NA)
  }  
dev.off()


## scale collapse
pdf('./figs/scale_collapse.pdf', width = 7 * 2, height = 7)
  ddr$ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)
  xlab = 'log(N/S) / log(Nbar/Sbar)'
  par(mfrow=c(1, 2))
  plot(b0 ~ ratio , data=ddr, xlab=xlab, ylab='b0', type='n')
  for (g in seq_along(grains)) 
    points(b0 ~ ratio, data=ddr, subset=grains == grains[g], col=g)
  legend('topright', legend=c('Grains', grains), col=c(NA, 1:length(grains)),
          pch=1, bty='n', cex = 2, lwd=3, lty=NA)
  plot(b1 ~ ratio , data=ddr, xlab=xlab, ylab='b1', type='n')
  for (g in seq_along(grains)) 
    points(b1 ~ ratio, data=ddr, subset=grains == grains[g], col=g)
dev.off()

## scale collapse lines
pdf('./figs/scale_collapse_lines.pdf', width = 7 * 2, height = 7*2)
  ddr$ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)
  xlab = 'log(N/S) / log(Nbar/Sbar)'
  par(mfrow=c(2, 2))
  plot(b0 ~ ratio , data=ddr, main='S fixed, N varies', xlab=xlab, ylab='b0', type='n')
  for (g in seq_along(grains)) {
    for (s in seq_along(S)) 
      lines(b0 ~ ratio, data=ddr, subset=grains == grains[g] & S == S[s], col=g)
  }
  plot(b0 ~ ratio , data=ddr, main='S varies, N fixed',xlab=xlab, ylab='b0', type='n')
  for (g in seq_along(grains)) {
    for (n in seq_along(N)) 
      lines(b0 ~ ratio, data=ddr, subset=grains == grains[g] & N == N[n], col=g)
  }
  legend('topright', legend=c('Grains', grains), col=c(NA, 1:length(grains)),
          pch=1, bty='n', cex = 2, lwd=3, lty=NA)
  plot(b1 ~ ratio , data=ddr, main='S fixed, N varies', xlab=xlab, ylab='b1', type='n')
  for (g in seq_along(grains)) {
    for (s in seq_along(S))
      lines(b1 ~ ratio, data=ddr, subset=grains == grains[g] & S == S[s], col=g)
  }  
  plot(b1 ~ ratio , data=ddr,  main='S varies, N fixed',xlab=xlab, ylab='b1', type='n')
  for (g in seq_along(grains)) {
    for (n in seq_along(N))
      lines(b1 ~ ratio, data=ddr, subset=grains == grains[g] & N == N[n], col=g)
  }  
dev.off()

## consider 3 of the grains
pdf('./figs/scale_collapse_3grains_pwr_wtr.pdf', width = 7 * 2,
    height = 7 * 2)
  ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)
  plot(ratio, pwrStats[cof, meth, g, , ], xlab=xlab, ylab=cof, type='n')
  xlab = 'log(N/S) / log(Nbar/Sbar)'
  lwd =2
  par(mfrow=c(3, 2))
  for (g in c(1, 3, 5)) {
    plot(b0 ~ ratio, data=ddr, subset = grains == grains[g], type='n',
         xlab=xlab, ylab='b0', 
         main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))
    for(s in seq_along(S))
      lines(b0 ~ ratio, data=ddr, subset = grains == grains[g] & S == S[s],
            col='palevioletred', lwd=lwd)
    plot(b0 ~ ratio, data=ddr, subset = grains == grains[g], type='n',
         xlab=xlab, ylab='b0', 
         main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))  
    for(n in seq_along(N))
      lines(b0 ~ ratio, data=ddr, subset = grains == grains[g] & N == N[n],
            col='dodgerblue', lwd=lwd)
  }
  for (g in c(1, 3, 5)) {
    plot(b1 ~ ratio, data=ddr, subset = grains == grains[g], type='n',
         xlab=xlab, ylab='b1', 
         main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))
    for(s in seq_along(S))
      lines(b1 ~ ratio, data=ddr, subset = grains == grains[g] & S == S[s],
            col='palevioletred', lwd=lwd)
    plot(b1 ~ ratio, data=ddr, subset = grains == grains[g], type='n',
         xlab=xlab, ylab='b1', 
         main = paste('Grain is', grains[g], 'of', 2^12, sep=' '))  
    for(n in seq_along(N))
      lines(b1 ~ ratio, data=ddr, subset = grains == grains[g] & N == N[n],
            col='dodgerblue', lwd=lwd)
  }  
dev.off()

## scale collapse for presentation
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  ddr$ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)
  xlab = 'log(N/S) / log(Nbar/Sbar)'
  par(mfrow=c(1, 1))
  plot(10^b0 ~ ratio , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(0, 1))
  axis(side=1, cex.axis=1.5, lwd=3)
  axis(side=2, cex.axis=1.5, lwd=3)
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio, 10^tmp$b0), col=col[g], lwd=4)
  }  
##
  plot(b1 ~ ratio , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F)
  axis(side=1, cex.axis=1.5, lwd=3)
  axis(side=2, cex.axis=1.5, lwd=3)  
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio, tmp$b1, f= 1/3), col=col[g], lwd=4)
  }  

mk_legend('center', legend=grains, col=col,
          lty=1, bty='n', cex = 2, lwd=8)


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


