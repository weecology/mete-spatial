setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

## Figure 2 - example Simulated DDR with simulation results---------------------
windows(width=7 * 3, height=7)
par(mfrow=c(1, 3))

axislwd = 4
linelwd = 3

## Example Simulated Distance Decay
## the simulated empirical data is used here instead of an example from the
## parameter space exploration b/c the Sherman simulated points have been binned
## in an attractive way. The simulated parameter space results are raw results
## witout binning
load('./simulated_empirical_results.Rdata')

tmp = list(simSorAbuLogSer$bormann_C200_B12_grid)
col = colorRampPalette(c('dodgerblue', 'red'))(5)
range(tmp[[1]]$Avg)
grains = unique(tmp[[1]]$Comm)

tmp[[1]]$Dist = tmp[[1]]$Dist / sqrt(grains[1])

plot(1:10, type='n', xlab='', ylab='', xlim=c(1, 50), 
     ylim = c(0.02, 0.45) , frame.plot=F, axes=F, log='xy')
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd)
axis(side=2, cex.axis=1.75, lwd=axislwd, at=c(0.02, 0.05, 0.10, 0.20, 0.40))
plotEmpir(tmp, 'average', log='xy', title=F, 
          quants=F, col=col, lwd=1, add=TRUE, type='p', cex=2, pch=19)
for (g in seq_along(grains)) {
  mod = lm(log(Avg) ~ log(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  lines(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]], exp(predict(mod)),
        col=col[g], lwd=linelwd)  
}

## these results were calculated in the script spat_param_space.R
ddr = read.csv('./sorensen/param_ddr_wtr_pwr_stats.csv')

## scale collapse for presentation
ddr$ratio = (log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)) * log(grains)
plot(10^b0 ~ ratio, data=ddr, subset=grains == 1, col='dodgerblue', pch=19,
     cex = 1.5, xlab='', ylab='', frame.plot=F, axes=F, xlim=c(0, 200))
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd, at = c(0, 50, 100, 150, 200))
axis(side=2, cex.axis=1.75, lwd=axislwd)
##
plot(b1 ~ ratio , data=ddr, subset=grains == 1, col='dodgerblue', pch=19,
     cex = 1.5, xlab='', ylab='', frame.plot=F, axes=F, xlim=c(0, 200))
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd, at = c(0, 50, 100, 150, 200))
axis(side=2, cex.axis=1.75, lwd=axislwd) 


## Figure 3: 6 panel SAR & DDR graphic------------------------------------------
## A) example SAR, B) SAR METE residuals, C) SAR RP residuals
## D) exampld DDR, E) DDR METE residuals, F) DDR RP residuals
##
## panel A) example SAR pattern
## load sar data and compute residuals
source('./scripts/spat_sar_load_and_avg_data.R')
## set up graphic parameters 
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
        "#FF8000")
windows(width=7*3, height=7*2)
par(mfrow=c(2,3))

## plot sar results
site = 'bigoak'
i = match(site, names(meteEmpirSAD))
  ## log-log
    plot(sr_iter ~ area, data=meteEmpirSAD[[i]], axes=F,
         ylim=range(c(60, meteEmpirSAD[[i]]$sr_iter, empir[[i]]$richness)),
         xlim=range(c(1, meteEmpirSAD[[i]]$area, empir[[i]]$area)), log='xy',
         type='n', ylab='', xlab='')
    addAxis1()
    addAxis2(at = 2^(-1:5))
    ## meteEmpirSAD CI
    dat = meteAvgEmpirSAD[[match(names(meteEmpirSAD)[i], names(meteAvgEmpirSAD))]]
    lines(sr.avg ~ grains, data=dat, lwd=3, lty=1, col=1)
    ## RP CI
    dat = srExp[[match(names(meteEmpirSAD)[i], names(srExp))]]
    lines(S_binom ~ grains, data=dat, col='grey', lty=1, lwd=3)
    ## analytical meteEmpirSAD    
#    lines(sr_noniter ~ area, data=meteEmpirSAD[[i]], col='dodgerblue', lwd=3)
    ## data
    lines(richness ~ area, data = empir[[i]], type='p', pch=19, lwd=3, cex=2)
    legend('bottomright', c('Observed', 'METE', 'RP'), pch=c(19, NA, NA), lwd=c(NA,5,5),
           col=c(1, 1, 'grey'), cex=2, bty='n')

## panels C & D: empirical SAR residuals
  sites = unique(sar_res$site)
  plot(empirsad_avg ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6))
  addAxis1(at=10 ^ seq(-1, 5, 2))
  addAxis2()
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_avg ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=habcol[habindex])
  }
#  legend('bottomright', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')
  ##
  plot(empirsad_avg ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6))
  addAxis1(at=10 ^ seq(-1, 5, 2))
  addAxis2()
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_rp ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=habcol[habindex], lty=1)
  }

## panel D) empirical DDR pattern at a single scale
load('./sorensen/empirSorAbu.Rdata')
load('simulated_empirical_results.Rdata')
## start with a good fit for slope by METE
tmp = empirSorAbu['bigoak']
obs = empirSorAbu$'bigoak'
exp = simSorAbuLogSer$'bigoak'
grains = unique(obs$Comm)

plot(Metric.avg ~ Dist, data=obs,
     ylim=c(0.02, .5), xlim=c(2,128), type='n', frame.plot=F, axes=F,
     xlab='', ylab='', log='xy')
addAxis1(at = 2^(1:7))
addAxis2(at= 0.02 * 2^(0:4))
g = 2
  ## add mete  
  true = exp$Comm == grains[g]
  tmpexp = exp[true,]
#  addCI('Dist', 'Avg.lo', 'Avg.hi', col='grey', data='tmpexp')
  lines(Avg ~ Dist, data=tmpexp, col=1, lty=1, lwd=3)
  ## RP
  lines(Exp.avg ~ Dist, data=obs, subset=Comm == grains[g], lty=1,
        col='grey', lwd=3)  
  ## data
  lines(Metric.avg ~ Dist, data=obs, subset=Comm == grains[g],
        lty=1, lwd=3, col=1, type='p', pch=19, cex=2)


#mk_legend('center', c('Observed', 'METE', 'Random Placement'),
#          col = c(1, lightblue, purple), lwd=3, lty = c(1,2,2),
#          cex=2, bty='n')

## panels E & F) emprical DDR residuals
resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)

dat = resSorAbuFixed
dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
sites = unique(dat$site)

  for(j in 1:2){
    if(j == 1)
      main = 'METE'
    else
      main = 'RP'
    plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim=c(-.05,.4), xlim=c(.5,512),
         xlab='', ylab='', axes=F, frame.plot=F)
    addAxis1(at=2^seq(-1, 9, 2))
    addAxis2()
    for(i in seq_along(sites)) {
      tmp = subset(dat, site == sites[i])
      grains = unique(tmp$Comm)
      habindex = match(habitat[match(sites[i], shrtnm)], hab)
      if(j == 1) {
          lines(lowess(tmp$Dist, tmp$avg.res), col = habcol[habindex], lwd=3)
      }
      else {
          lines(lowess(tmp$Dist, tmp$exp.res), col = habcol[habindex], lwd=3)  
      }  
    }  
  }  
#mk_legend('center', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')

## METE DDR scale colapse with data---------------------------------------------
empirStatsSorAbuAvg = getStats(empirSorAbu, 'average')

fileNames = dir('./sar')
empirFiles = grep('empir_sar.csv', fileNames)
empir = vector('list', length(empirFiles))
names(empir) = sub('_empir_sar.csv', '', fileNames[empirFiles])
for (i in seq_along(empirFiles)) {
  empir[[i]] = read.csv(paste('./sar/', fileNames[empirFiles[i]], sep=''))
  empir[[i]]$area = round(empir[[i]]$area, 2)
}

## convert to flat file
stats = empirStatsSorAbuAvg
mod = 'pwr'
meth = 'wtr'
rm(dd_stats)
for (i in seq_along(stats)) {
  site = names(stats)[i]
  index = match(site, names(empir))
  grains = as.numeric(dimnames(stats[[i]])[[4]])
  S = max(empir[[index]]$richness)
  N = max(empir[[index]]$indiv)
  for (g in seq_along(grains)) {
    empir_tmp = empir[[index]][empir[[index]]$area == grains[g], ]
    b0 = stats[[i]][mod, 'b0', meth, g]
    b1 = stats[[i]][mod, 'b1', meth, g]
    navg = empir_tmp$indiv
    savg = empir_tmp$richness
    ratio = log(N / S) / log(navg / savg)
    if (exists('dd_stats'))
      dd_stats = rbind(dd_stats, 
                       data.frame(site, area=grains[g], S, N, savg, navg, b0, b1, ratio))
    else
      dd_stats = data.frame(site, area=grains[g], S, N, savg, navg, b0, b1, ratio)

  }
}

## bring in scale collapse information from the parameter space analysis

sites = unique(dd_stats$site)

par(mfrow=c(1,2))
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  sar$ratio = log(sar$N / sar$S) / log(sar$ind.avg / sar$sr.avg)
  xlab = 'log(N/S) / log(Nbar/Sbar)'
  plot(10^b0 ~ ratio , data=sar, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(0, 1), xlim=c(0,80))
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8)
  axis(side=2, cex.axis=1.75, lwd=8)
  for (g in seq_along(grains)) { 
    tmp = subset(sar, grains == grains[g])
    lo = lowess(tmp$ratio, 10^tmp$b0)
    true = lo$y < 1 & lo$y >0
    lines(lo$x[true], lo$y[true], col='grey',lty=3, lwd=4)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(10^b0 ~ ratio, data=dd_stats, subset= site == sites[i], 
           col=habcol[habindex], lwd=4)
  }
##
  plot(b1 ~ ratio , data=sar, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, xlim=c(0, 80), ylim=c(-.65, 0))
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8)
  axis(side=2, cex.axis=1.75, lwd=8) 
  for (g in seq_along(grains)) {
    tmp = subset(sar, grains == grains[g])
    lines(lowess(tmp$ratio, tmp$b1, f= 1/3), col='grey', lty=3, lwd=4)
  }
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(b1 ~ ratio, data=dd_stats, subset= site == sites[i], 
           col=habcol[habindex], lwd=4)
  }



## Supplemental Figure 1--------------------------------------------------------
## r2 of model fits to simulated results
#load('./sorensen/simSorAbuAvg.Rdata')
#stats = getSimStats(simSorAbuAvg, S, N)
windows(width= 7 * 1, height=7)
par(mfrow=c(1,1))

meth='ols'
dpwr = density(stats['pwr', 'r2', meth, , , ], na.rm = TRUE)
dexp = density(stats['exp', 'r2', meth, , , ], na.rm = TRUE)
dpwr$x = c(.7, dpwr$x)
dpwr$y = c(0, dpwr$y)
## drop x values less than 0.7
dpwr$y = dpwr$y[dpwr$x >= .7]
dpwr$x = dpwr$x[dpwr$x >= .7]
dexp$y = dexp$y[dexp$x >= .7]
dexp$x = dexp$x[dexp$x >= .7]

xlims = range(c(dpwr$x, dexp$x, 1))
ylims = range(c(dpwr$y, dexp$y))

plot(dpwr$x, dpwr$y, type='l', lty=3, lwd=linelwd, xlim=round(xlims,1), ylim=ylims, col='black',
     xlab='', ylab='', frame.plot=F, axes=F)
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd, at=c(.7, .8, .9, 1))
#axis(side=2, cex.axis=1.75, lwd=axislwd)
lines(dexp$x[dexp$x <=1.01], dexp$y[dexp$x <=1.01], lty=3, lwd=linelwd, col='grey')
##
meth='wtr'
dpwr = density(stats['pwr', 'r2', meth, , , ], na.rm = TRUE)
dexp = density(stats['exp', 'r2', meth, , , ], na.rm = TRUE)
dpwr$x = c(min(dexp$x), dpwr$x)
dpwr$y = c(0, dpwr$y)
xlims = range(c(dpwr$x, dexp$x, 1))
ylims = range(c(dpwr$y, dexp$y))

lines(dpwr$x, dpwr$y, lwd=linelwd, col='black')
lines(dexp$x[dexp$x <=1.01], dexp$y[dexp$x <=1.01], lwd=linelwd, col='grey')

mk_legend('center', c('Exponential, OLS', 'Exponential, WLS', 'Power, OLS', 'Power, WLS'),
          lty=c(3, 1, 3, 1), lwd=6, bty='n', cex=2, col=c('grey', 'grey', 'black','black'))

## Suplemental scale collapse figure with all grains----------------------------
## these results were calculated in the script spat_param_space.R
ddr = read.csv('./sorensen/param_ddr_wtr_pwr_stats.csv')

## scale collapse for presentation
par(mfrow=c(1,2))
ddr$ratio = (log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)) * log(grains)
xlab = 'log(N/S) / log(Nbar/Sbar)'
plot(10^b0 ~ ratio , data=ddr, xlab='', ylab='', type='n',
     frame.plot=F, axes=F, ylim=c(0, 1), xlim=range(ddr$ratio))
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd)
axis(side=2, cex.axis=1.75, lwd=axislwd)
for (g in seq_along(grains)) { 
  tmp = subset(ddr, grains == grains[g])
  lo = lowess(tmp$ratio, 10^tmp$b0)
#  true = lo$y < 1 & lo$y >0
#  lines(lo$x[true], lo$y[true], col=col[g], lwd=linelwd)
  lines(lo$x, lo$y, col=col[g], lwd=linelwd)
}  
##
plot(b1 ~ ratio , data=ddr, xlab='', ylab='', type='n',
     frame.plot=F, axes=F, xlim=range(ddr$ratio))
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd)
axis(side=2, cex.axis=1.75, lwd=axislwd) 
for (g in seq_along(grains)) {
  tmp = subset(ddr, grains == grains[g])
  lines(lowess(tmp$ratio, tmp$b1, f= 1/3), col=col[g], lwd=linelwd)
}  
##
mk_legend('center', legend=grains, col=col,
          lty=1, bty='n', cex = 2, lwd=8)




