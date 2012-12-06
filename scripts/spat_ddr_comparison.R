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
par(mfrow=c(1,2))
for(k in 1:2){
  for(j in 1:2){
    if(j == 1)
      main = 'METE'
    else
      main = 'RP'
    plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim=c(-.05,.4), main=main)
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
    brks = seq(0, 1, length.out=21)
    par(mfrow=c(1,2))
    hist(r2exp, xlim=c(0,1), breaks=brks, col='red')
    hist(r2pwr, xlim=c(0,1), breaks=brks, col='blue', add=TRUE)
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
