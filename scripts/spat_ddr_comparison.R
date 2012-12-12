## Purpose: to compare the fit of the simulated ddr to the empirical ddr

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

## load the simulated and empirical DDR results---------------------------------
load('./simulated_empirical_results.Rdata')
load('./sorensen/empirSorBin.Rdata') 
load('./sorensen/empirSorAbu.Rdata') 
load('./varWithin/empirVarBin.Rdata')
load('./varWithin/empirVarAbu.Rdata')

## compute residuals between observed and expected-------------------------------
resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)
resSorAbuLogSer = get_ddr_resid(empirSorAbu, simSorAbuLogSer)
resVarAbuFixed = get_ddr_resid(empirVarAbu, simVarAbuFixed)
resVarAbuLogSer = get_ddr_resid(empirVarAbu, simVarAbuLogSer)
##
resSorBinFixed = get_ddr_resid(empirSorBin, simSorBinFixed)
resSorBinLogSer = get_ddr_resid(empirSorBin, simSorBinLogSer)
resVarBinFixed = get_ddr_resid(empirVarBin, simVarBinFixed)
resVarBinLogSer = get_ddr_resid(empirVarBin, simVarBinLogSer)

## examine empirical fit for a single dataset------------------------------------
tmp = empirSorAbu['bigoak']
plotEmpir(tmp, 'average', log='xy')
obs = empirSorAbu$'bigoak'
exp = simSorAbuFixed$'bigoak'
grains = unique(obs$Comm)

purple = rgb(112, 48, 160, maxColorValue=160)
lightblue = "#1AB2FF"

plot(Metric.avg ~ Dist, data=obs,
     ylim=c(0.02, .5), xlim=c(2,128), type='n', frame.plot=F, axes=F,
     xlab='', ylab='', log='xy')
axis(side=1, cex.axis=1.75, lwd=8, padj=.5,
     at = 2^(1:7))
axis(side=2, cex.axis=1.75, lwd=8,
     at= 0.02 * 2^(0:4))
for (g in 2) {
  ## add mete  
  true = exp$Comm == grains[g]
  tmpexp = exp[true,]
  addCI('Dist', 'Avg.lo', 'Avg.hi', data='tmpexp', col='grey')
  lines(Avg ~ Dist, data=tmpexp, col=lightblue, lty=2, lwd=4)
  ## RP
  true = obs$Comm == grains[g]
  tmpobs = obs[true, ]
  addCI('Dist', 'Exp.25', 'Exp.75', data='tmpobs', col='grey')
  lines(Exp.avg ~ Dist, data=obs, subset=Comm == grains[g], lty=2,
        col=purple, lwd=4)  
  ## data
  lines(Metric.avg ~ Dist, data=obs, subset=Comm == grains[g],
        lty=1, lwd=4, col=1, type='o', pch=19, cex=1.25)  
}

mk_legend('center', c('Observed', 'METE', 'Random Placement'),
          col = c(1, lightblue, purple), lwd=3, lty = c(1,2,2),
          cex=2, bty='n')

## Examine residuals------------------------------------------------------------

plot(med.res ~ Dist, data=resSorAbuLogSer, log='x')

par(mfrow=c(1,2))
plot(resSorAbuFixed$avg.res[resSorAbuFixed$site != 'cross'],
     resSorAbuLogSer$avg.res)
abline(a=0, b=1)
sum(resSorAbuFixed$avg.res[resSorAbuFixed$site != 'cross'])
sum(resSorAbuLogSer$avg.res)
diff = resSorAbuFixed$avg.res[resSorAbuFixed$site != 'cross'] - 
       resSorAbuLogSer$avg.res
plot(resSorAbuLogSer$Dist, diff, log='x')
abline(a=0,b=0)

## choose which simulated result is of most interest
dat = resSorAbuFixed
dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
sar_data = read.csv('./sar/empir_sars.csv')
sar_data$area = round(sar_data$area, 2)
dat = merge(dat, sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)
## subset so that has at least 20 individuals
dat = subset(dat, indiv >= 20)
## normalize by S
#dat[, 14:16] = dat[ , 14:16]  / dat$richness

sites = unique(dat$site)

avgresid = aggregate(avg.res ~ site, data=dat, function(x) sum(x^2) / length(x))
avgresid[order(avgresid[ , 2]), ]
expresid = aggregate(exp.res ~ site, data=dat, function(x) sum(x^2) / length(x))
expresid[order(expresid[ , 2]), ]

plot(avgresid[,2], expresid[,2])
abline(a=0, b=1)

shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
         "#FF8000")

## residuals vs distance each grain has its own line
par(mfrow=c(1,1))
plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim=c(-.1,.6))
abline(h=0, lty=2, lwd=2)
for(i in seq_along(sites)) {
  tmp = subset(dat, site == sites[i])
  grains = unique(tmp$Comm)
  col = colorRampPalette(c('dodgerblue', 'red'))(length(grains))
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  for(g in seq_along(grains)) {
    lines(avg.res ~ Dist, data=tmp, subset=Comm==grains[g],
          col = habcol[habindex], lwd=2)
    lines(exp.res ~ Dist, data=tmp, subset=Comm==grains[g],
          col = habcol[habindex], lwd=2, lty=2)
  }  
}
mod = lm(avg.res ~ log(Dist), data=dat)
abline(mod)
summary(mod)

## residuals vs distance summarized per dataset via lowess lines
pdf('./figs/sor_abu_fixed_residuals.pdf', width = 7 *2 , height=7)
ylims = c(-.05,.6) #c(-.05, .3)
par(mfrow=c(1,2))
for(k in 1:2){
  for(j in 1:2){
    if(j == 1)
      main = 'METE'
    else
      main = 'RP'
    plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim=ylims , main=main)
    abline(h=0, lty=2, lwd=2)
    for(i in seq_along(sites)) {
      tmp = subset(dat, site == sites[i])
      grains = unique(tmp$Comm)
      col = colorRampPalette(c('dodgerblue', 'red'))(length(grains))
      habindex = match(habitat[match(sites[i], shrtnm)], hab)
      if(j == 1) {
        if(k == 1)
          points(tmp$Dist, tmp$avg.res, col = habcol[habindex], pch=19)
        else
          lines(lowess(tmp$Dist, tmp$avg.res), col = habcol[habindex], lwd=2)
      }
      else {
        if(k ==1)
          points(tmp$Dist, tmp$exp.res, col = habcol[habindex], pch=19)
        else
          lines(lowess(tmp$Dist, tmp$exp.res), col = habcol[habindex], lwd=2, lty=2)  
      }  
    }  
  }
}  
mk_legend('center', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')
dev.off()

## residuals vs area
par(mfrow=c(1,1))
plot(avg.res ~ area, data=dat, type='n', log='x', axes=F, frame.plot=F,
     xlim = c(0.1, 1e5), ylim=c(0, .6), xlab='', ylab='')
axis(side=1, cex.axis=1.75, padj=.5, lwd=8,
     at=10 ^ (-1:5))
axis(side=2, cex.axis=1.75, lwd=8)
for(i in seq_along(sites)) {
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  tmp = subset(dat, site == sites[i])
  lines(lowess(tmp$area, tmp$avg.res),
         col=habcol[habindex], pch=1, lwd=4)
  lines(lowess(tmp$area, tmp$exp.res), lty=2,
         col=habcol[habindex], pch=1, lwd=4)  
}
mk_legend('center', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')


## avg across the distances 
resMETE  = aggregate(dat$avg.res, 
                 by=list(site=dat$site, area=dat$area),
                 FUN=function(x) sum(x^2)/length(x))
names(resMETE ) = c('site','area','res')

resRP  = aggregate(dat$exp.res, 
                 by=list(site=dat$site, area=dat$area),
                 FUN=function(x) sum(x)/length(x))
names(resRP) = c('site','area','res')

plot(res ~ area, data=resMETE , log='x', type='n',
     ylim=c(-.4, .3))
for(i in seq_along(sites)) {
  tmp = subset(resMETE , )
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  lines(res ~ area, data=resMETE, subset=site == sites[i],
        col=habcol[habindex], pch=1, lwd=2)
  lines(res ~ area, data=resRP, subset=site == sites[i],
        col=habcol[habindex], pch=1, lwd=2, lty=2)
}
abline(h=0)


mk_legend('center', hab, col=col, pch=19)
  
## compute summary statics for empricial and simualted results-------------------
empirStatsSorAbuMed = getStats(empirSorAbu, 'median')
empirStatsSorAbuAvg = getStats(empirSorAbu, 'average')
  
simStatsSorAbuFixedMed = getStats(simSorAbuFixed, 'median')
simStatsSorAbuFixedAvg = getStats(simSorAbuFixed, 'average')

simStatsSorAbuLogSerMed = getStats(simSorAbuLogSer, 'median')
simStatsSorAbuLogSerAvg = getStats(simSorAbuLogSer, 'average')

##
empirStatsVarAbuMed = getStats(empirVarAbu, 'median')
empirStatsVarAbuAvg = getStats(empirVarAbu, 'average')

simStatsVarAbuFixedMed = getStats(simVarAbuFixed, 'median')
simStatsVarAbuFixedAvg = getStats(simVarAbuFixed, 'average')

simStatsVarAbuLogSerMed = getStats(simVarAbuLogSer, 'median')
simStatsVarAbuLogSerAvg = getStats(simVarAbuLogSer, 'average')

##-----------------------------------------------------------------------------

empirStatsSorBinMed = getStats(empirSorBin, 'median')
empirStatsSorBinAvg = getStats(empirSorBin, 'average')

simStatsSorBinFixedMed = getStats(simSorBinFixed, 'median')
simStatsSorBinFixedAvg = getStats(simSorBinFixed, 'average')

simStatsSorBinLogSerMed = getStats(simSorBinLogSer, 'median')
simStatsSorBinLogSerAvg = getStats(simSorBinLogSer, 'average')

##
empirStatsVarBinMed = getStats(empirVarBin, 'median')
empirStatsVarBinAvg = getStats(empirVarBin, 'average')

simStatsVarBinFixedMed = getStats(simVarBinFixed, 'median')
simStatsVarBinFixedAvg = getStats(simVarBinFixed, 'average')

simStatsVarBinLogSerMed = getStats(simVarBinLogSer, 'median')
simStatsVarBinLogSerAvg = getStats(simVarBinLogSer, 'average')

## generate pdfs of r2 values, exponential or pwr function better
pl_r2 = function() {
  for( i in 1:2) {
    r2pwr = unlist(sapply(stats[[i]], function(x) x['pwr', 'r2', 'wtr',]))
    r2exp = unlist(sapply(stats[[i]], function(x) x['exp', 'r2', 'wtr',]))
    dpwr = density(r2pwr, na.rm=T)
    dexp = density(r2exp, na.rm=T)
    par(mfrow=c(1,2))
    plot(dpwr, col='blue',
         xlim=range(c(dpwr$x, dexp$x)), ylim=range(c(dpwr$y, dexp$y)))
    lines(dexp, col='red')
    legend('topleft', c('pwr','exp'), col=c('blue','red'), lty=1)
    plot(r2pwr, r2exp)
    abline(a=0, b=1)
  }    
}

stats = list(empirStatsSorAbuAvg,  simStatsSorAbuFixedAvg)
pdf('./figs/exp_vs_pwr_avg_sor_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsSorBinAvg,  simStatsSorBinFixedAvg)
pdf('./figs/exp_vs_pwr_avg_sor_bin_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsSorAbuMed,  simStatsSorFixedMed)
pdf('./figs/exp_vs_pwr_med_sor_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsVarAbuAvg, simStatsVarAbuFixedAvg, simStatsVarAbuLogSerAvg)
pdf('./figs/exp_vs_pwr_avg_var_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsVarAbuMed, simStatsVarAbuFixedMed, simStatsVarAbuLogSerMed)
pdf('./figs/exp_vs_pwr_med_var_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  

## plot summary statistics------------------------------------------------------

pl_coef = function(stats1, stats2, stats3=NULL, xlim=NULL, ylim=NULL) {
  if (is.null(xlim))
    xlim = range(unlist(sapply(stats1, function(x) range(x[mod, cof, meth, ], na.rm=T))))
  if (is.null(ylim))
    ylim = range(unlist(sapply(stats2, function(x) range(x[mod, cof, meth, ], na.rm=T))))
  plot(stats1$bci[mod, cof, meth, ], stats2$bci[mod, cof, meth, ], 
       type='n', xlim = xlim, ylim= ylim, main=paste(mod, meth, sep=', '),
       xlab=paste('Empir', cof, sep=' '), ylab=paste('METE', cof, sep=' '))
  for (i in seq_along(stats1)) {
    site = names(stats1)[i]
    index2 = which(grepl(site, names(stats2)))
    if (length(index2) > 0)
      points(stats1[[i]][mod, cof, meth, ], 
             stats2[[index2]][mod, cof, meth, ], col=col[i], pch=19, cex=1.25, type='o')
    if (!is.null(stats3)) {
      index3 = which(grepl(site, names(stats3)))
      if (length(index3) > 0)
        points(stats1[[i]][mod, cof, meth, ], 
               stats3[[index3]][mod, cof, meth, ], col=col[i])
    }  
  }
  abline(a=0, b=1)
}

pdf('./figs/coef_one_to_one_sor_abu_bin_avg_wtr.pdf', width=7 * 2, height=7*2)
  stats1 = empirStatsSorAbuAvg
  stats2 = simStatsSorAbuFixedAvg
  stats1 = stats1[-grep('serp', names(stats1))]
  stats2 = stats2[-grep('serp', names(stats2))]

  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  #col = terrain.colors(length(stats1) + 2)
  col = colorRampPalette(c("blue", "red", "green","orange"))(length(stats1))
  pl_coef(stats1, stats2)
  legend('topleft', names(stats1), cex=2, col=col, pch=19, bty='n')
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)

##------------------------------
  stats1 = empirStatsSorBinAvg
  stats2 = simStatsSorBinFixedAvg
  stats1 = stats1[-grep('serp', names(stats1))]
  stats2 = stats2[-grep('serp', names(stats2))]

  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  legend('topleft', names(stats1), cex=2, col=col, pch=19, bty='n')

  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
dev.off()

##-----------------------------------------------------------------------------
pdf('./figs/coef_one_to_one_var_abu_bin_avg_wtr.pdf', width=7 * 2, height=7*1.5)
  stats1 = empirStatsVarAbuAvg 
  stats2 = simStatsVarAbuFixedAvg
  stats3 = simStatsVarAbuLogSerAvg
  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
##------------------------------
  stats1 = empirStatsVarBinAvg 
  stats2 = simStatsVarBinFixedAvg
  stats3 = simStatsVarBinLogSerAvg
  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
dev.off()  

##-----------------------------------------------------------------------------
## examine scale collapse from simualted parameter space with empirical data----
ddr = read.csv('./sorensen/param_ddr_wtr_pwr_stats.csv')
## read in stats for the empirical datasets
load('./sorensen/empirSorAbu.Rdata')
empirStatsSorAbuAvg = getStats(empirSorAbu, 'average')

fileNames = dir('./sar')
empirFiles = grep('empir_sar.csv', fileNames)
empir = vector('list', length(empirFiles))
names(empir) = sub('_empir_sar.csv', '', fileNames[empirFiles])
for (i in seq_along(empirFiles)) {
  empir[[i]] = read.csv(paste('./sar/', fileNames[empirFiles[i]], sep=''))
  ## add rounded areas for lookup matching purposes
  empir[[i]]$area_r = round(empir[[i]]$area, 2)
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
  ## get areas unrounded
  areas = empir[[index]]$area[ empir[[index]]$area_r %in% grains]
  S0 = max(empir[[index]]$richness)
  N0 = max(empir[[index]]$indiv)
  A0 = max(empir[[index]]$area)
  for (g in seq_along(grains)) {
    empir_tmp = empir[[index]][empir[[index]]$area_r == grains[g], ]
    b0 = stats[[i]][mod, 'b0', meth, g]
    b1 = stats[[i]][mod, 'b1', meth, g]
    navg = empir_tmp$indiv
    savg = empir_tmp$richness
    A = areas[g]
    ratio_b0 = log(N0 / S0) / log(navg / savg)
    ratio_b1 = log(N0 / S0) / log(navg / savg) / log(A0 / A)
    if (exists('dd_stats'))
      dd_stats = rbind(dd_stats, 
                       data.frame(site, area=grains[g], S, N, savg, navg, b0, b1, ratio_b0, ratio_b1))
    else
      dd_stats = data.frame(site, area=grains[g], S, N, savg, navg, b0, b1, ratio_b0, ratio_b1)
  }
}

## bring in scale collapse information from the parameter space analysis

sites = unique(dd_stats$site)

## set up graphic parameters 
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
        "#FF8000")

#pdf('./figs/mete_ddr_scale_collapse_with_empir_data.pdf', width= 7 * 2, height=7 * 1)
#windows(width = 7 * 2, height=7 * 1)

ddr$ratio_b0 = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)
ddr$ratio_b1 = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg) / log(4096 / ddr$grains)

par(mfrow=c(1,2))
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  plot(10^b0 ~ ratio_b0 , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, log='', ylim=c(0, 2))
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(bar(italic(N))/bar(italic(S)))), padj=2.25)
  addylab('Decay Rate')
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    points(tmp$ratio_b0, 10^tmp$b0, col='grey', lwd=1, pch=19)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    points(10^b0 ~ ratio_b0, data=dd_stats, subset= site == sites[i] & navg >= 20, 
           col=habcol[habindex], lwd=2, lty=2)
  }
##
  plot(b1 ~ ratio_b1 , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(-.7, .1), log='')
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(bar(italic(N))/bar(italic(S))) / log(italic(A[0])/italic(A))), padj=2.25)
  addylab('Decay Rate')
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    points(tmp$ratio_b1, tmp$b1, col='grey', lwd=1, pch=19)
#    lines(lowess(tmp$ratio_b1, tmp$b1), col=col[g], lwd=4)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    points(b1 ~ ratio_b1, data=dd_stats, subset= site == sites[i] & navg >= 20, 
           col=habcol[habindex], lwd=2, lty=2)
  }
  legend('topright', c('METE Simulated', hab), cex=1.5,
         col=c('grey', habcol), lty=NA, lwd=3, pch=c(19, rep(1,5)), bty='n')

## Examine both possible ways of scaling the x-axis for the collapse-------------

pdf('./figs/mete_ddr_scale_collapse_both_scalings.pdf',
    width= 7 * 2, height=7 * 2)

par(mfrow=c(2,2))
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  plot(10^b0 ~ ratio_b0 , data=ddr, xlab='', ylab='', type='n',
       frame.plot=T, axes=F, log='', ylim=c(0, 1))
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(bar(italic(N))/bar(italic(S)))),
          padj=1.75)
  addylab('Decay Rate')
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio_b0, 10^tmp$b0), col=col[g], lwd=2)
  }  
##
  plot(b1 ~ ratio_b0 , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(-.7, .1), log='')
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(bar(italic(N))/bar(italic(S)))),
          padj=1.75)
  addylab('Decay Rate')
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio_b0, tmp$b1), col=col[g], lwd=2)
  }  
  legend('topright', legend=c('Grains',  unique(ddr$grains)), cex=2,
         col=c(NA, col), lwd=c(NA, rep(4, 5)), bty='n')
####
  plot(10^b0 ~ ratio_b1 , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, log='', ylim=c(0, 1))
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(bar(italic(N))/bar(italic(S))) / log(italic(A[0])/italic(A))),
          padj=1.75)
  addylab('Decay Rate')
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio_b1, 10^tmp$b0), col=col[g], lwd=2)
  }  
##
  plot(b1 ~ ratio_b1 , data=ddr, xlab='', ylab='', type='n',
       frame.plot=T, axes=F, ylim=c(-.7, .1), log='')
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(bar(italic(N))/bar(italic(S))) / log(italic(A[0])/italic(A))),
          padj=1.75)
  addylab('Decay Rate')
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio_b1, tmp$b1), col=col[g], lwd=2)
  }  
dev.off()

