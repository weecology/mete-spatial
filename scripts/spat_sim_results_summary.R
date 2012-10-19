## Purpose: to graphically summarize the simulated simical results
setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

files = get_file_names(path = './sorensen/', 
                      strings = list(shrtnames, 'C200', 'sherman3'), 
                      invert = c(F, F, T))

longnames = sub('sorensen_', '', files)
longnames = sub('_abu.Rdata', '', longnames)
longnames = unique(sub('_binary.Rdata', '', longnames))

simSorBin = getResults(longnames, 'sorensen', 'binary', TRUE)
simSorAbu = getResults(longnames, 'sorensen', 'abu', TRUE)

simSorBinAvg = avgSimResults(simSorBin, 'sorensen')
simSorAbuAvg = avgSimResults(simSorAbu, 'sorensen')

simVarBin = getResults(longnames, 'varWithin', 'binary', TRUE)
simVarAbu = getResults(longnames, 'varWithin', 'abu', TRUE)

simVarBinAvg = avgSimResults(simVarBin, 'varWithin')
simVarAbuAvg = avgSimResults(simVarAbu, 'varWithin')

## change names to better match empirical names
for (i in seq_along(longnames)) {
  newnames = c('Comm', 'Dist', 'N', 'Med.lo', 'Med', 'Med.hi', 'Avg.lo', 'Avg',
               'Avg.hi', 'Exp.lo', 'Exp', 'Exp.hi')
  names(simSorBinAvg[[i]]) = newnames
  names(simSorAbuAvg[[i]]) = newnames
  names(simVarBinAvg[[i]]) = newnames
  names(simVarAbuAvg[[i]]) = newnames
}

## export as files
#save(simSorBinAvg, simSorAbuAvg, file = './sorensen/simEmpirSorAvg.Rdata')
#save(simVarBinAvg, simVarAbuAvg, file = './varWithin/simEmpirVarAvg.Rdata')
load('./sorensen/simEmpirSorAvg.Rdata')
load('./varWithin/simEmpirVarAvg.Rdata')

## split results that were generated from log-series and those that were fixed abu
logser = grep('empirSAD', longnames, invert=TRUE)
fixed = grep('empirSAD', longnames, invert=FALSE)
simSorBinLogSer = simSorBinAvg[logser]
simSorBinFixed = simSorBinAvg[fixed]

simSorAbuLogSer = simSorAbuAvg[logser]
simSorAbuFixed = simSorAbuAvg[fixed]

simVarBinLogSer = simVarBinAvg[logser]
simVarBinFixed = simVarBinAvg[fixed]

simVarAbuLogSer = simVarAbuAvg[logser]
simVarAbuFixed = simVarAbuAvg[fixed]

## merge cocoli 1 & 2 and sherman 1 & 2
simSorBinLogSer = merge_drop(simSorBinLogSer)
simSorBinFixed = merge_drop(simSorBinFixed )

simSorAbuLogSer = merge_drop(simSorAbuLogSer)
simSorAbuFixed = merge_drop(simSorAbuFixed)

simVarBinLogSer = merge_drop(simVarBinLogSer)
simVarBinFixed = merge_drop(simVarBinFixed)

simVarAbuLogSer = merge_drop(simVarAbuLogSer)
simVarAbuFixed = merge_drop(simVarAbuFixed)

## update shrnames
shrtnames = shrtnames[-match(c('cocoli2','sherman2', 'sherman3'), shrtnames)] 

## export results
save(simSorAbuLogSer, simSorAbuFixed, simSorBinLogSer, simSorBinFixed,
     simVarAbuLogSer, simVarAbuFixed, simVarBinLogSer, simVarBinFixed,
     file = './simulated_empirical_results.Rdata')

## load the Results
load('./simulated_empirical_results.Rdata')
load('./sorensen/empirSorBin.Rdata') 
load('./sorensen/empirSorAbu.Rdata') 
load('./varWithin/empirVarBin.Rdata')
load('./varWithin/empirVarAbu.Rdata')

## drop ferp last grain
empirSorBin$ferp = empirSorBin$ferp[-nrow(empirSorBin$ferp),]
empirSorAbu$ferp = empirSorAbu$ferp[-nrow(empirSorAbu$ferp),]
empirVarBin$ferp = empirVarBin$ferp[-nrow(empirVarBin$ferp),]
empirVarAbu$ferp = empirVarAbu$ferp[-nrow(empirVarAbu$ferp),]



## examine a single simulated distance decay pattern
tmp = simSorAbuLogSer[14]
col = colorRampPalette(c('dodgerblue', 'red'))(5)
range(tmp[[1]]$Avg)
grains = unique(tmp[[1]]$Comm)

plot(1:10, type='n', xlab='', ylab='', xlim=c(0, 110), 
     ylim = c(0, 0.3) , frame.plot=F, axes=F, log='')
axis(side=1, cex.axis=1.25, lwd=2)
axis(side=2, cex.axis=1.25, lwd=2)
plotEmpir(tmp, 'average', log='xy', title=F, 
          quants=F, col=col, lwd=5, add=TRUE, type='p', pch=19)
for (g in seq_along(grains)) {
  mod = lm(log(Avg) ~ log(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  lines(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]], exp(predict(mod)),
        col=col[g], lwd=4)  
}

###
plot(1:10, type='n', xlab='', ylab='', xlim=c(1.5, 120), 
     ylim = c(0.006, 0.40) , frame.plot=F, axes=F, log='xy')
axis(side=1, cex.axis=1.75, padj=.5, lwd=8)
axis(side=2, cex.axis=1.75, lwd=8, at=c(0.01, 0.02, 0.05, 0.10, 0.20, 0.40))
plotEmpir(tmp, 'average', log='xy', title=F, 
          quants=F, colcol, lwd=5, add=TRUE, type='p', cex=1.5, pch=19)
for (g in seq_along(grains)) {
  mod = lm(log(Avg) ~ log(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  lines(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]], exp(predict(mod)),
        col=col[g], lwd=6)  
}
##
plot(log2(Avg) ~ log2(Dist), data=tmp[[1]], type='n', xlab='', ylab='',
     frame.plot=F, axes=F)
axis(side=1, cex.axis=1.25, lwd=2)
axis(side=2, cex.axis=1.25, lwd=2)
for (g in seq_along(grains)) {
  mod = lm(log2(Avg) ~ log2(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  points(log2(Avg) ~ log2(Dist), data=tmp[[1]], subset=Comm == grains[g],
         col=col[g], pch=19, lwd=5)
  lines(log2(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]]), predict(mod),
        col=col[g], lwd=4)  
}

##------------------------------------------------------------------------------
## examine empirical fit for a single dataset
tmp = empirSorAbu['bigoak']
plotEmpir(tmp, 'average', log='xy')
obs = empirSorAbu$'bigoak'
exp = simSorAbuLogSer$'bigoak'
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
  xvals = c(tmpexp$Dist, rev(tmpexp$Dist))
  yvals = c(tmpexp$Avg.lo, rev(tmpexp$Avg.hi))
  polygon(xvals, yvals, col="grey", border=NA)
  lines(Avg ~ Dist, data=tmpexp, col=lightblue, lty=2, lwd=4)
  ## RP
  lines(Exp ~ Dist, data=obs, subset=Comm == grains[g], lty=2,
        col=purple, lwd=4)  
  ## data
  lines(Metric.avg ~ Dist, data=obs, subset=Comm == grains[g],
        lty=1, lwd=4, col=1, type='o', pch=19, cex=1.25)  
}

mk_legend('center', c('Observed', 'METE', 'Random Placement'),
          col = c(1, lightblue, purple), lwd=3, lty = c(1,2,2),
          cex=2, bty='n')

##-------------------------------------------------------------------------------

## export pdf summaries of the empirical and simulated results
mk_summary = function() {
    par(mfrow=c(2,3))
    comp_results(site,
                 list(results1),
                 list("quants=F, expected=TRUE, sub='RandomPlacement',
                       metric='average', ylim=ylim, col=col, log=log_tran"))
    comp_results(site,
                 list(results1, results2),
                 list("quants=F, sub='Fixed'",
                      "quants=T, lwd=2, lty=2"),
                 metric='average', ylim=ylim, col=col, log=log_tran)
    comp_results(site,
                 list(results1, results3),
                 list("quants=F, sub='LogSer'",
                      "quants=T, lwd=2, lty=2"),
                 metric='average', ylim=ylim,col=col, log=log_tran)
    ###
    comp_results(site,
                 list(results1),
                 list("quants=F, ylim=ylim,
                       col=col, log=log_tran"))
    comp_results(site,
                 list(results1, results2),
                 list("quants=T, sub='AbuFixed'",
                      "quants=T, lwd=2, lty=2"),
                 metric='median', ylim=ylim, col=col, log=log_tran)
    comp_results(site,
                 list(results1, results3),
                 list("quants=T, sub='LogSer'",
                      "quants=T, lwd=2, lty=2"),
                 metric='median', ylim=ylim,col=col, log=log_tran)
}

results1 = empirSorAbu
results2 = simSorAbuFixed
results3 = simSorAbuLogSer
metric = 'sorensen'
dataType = 'abu'
log_tran='xy'
ylim=list(c(0.1,1))
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------

results1 = empirBinAbu
results2 = simSorBinFixed
results3 = simSorBinLogSer
metric = 'sorensen'
dataType = 'bin'
log_tran='xy'
ylim=list(c(0.1,1))
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------
results1 = empirVarAbu
results2 = simVarAbuFixed
results3 = simVarAbuLogSer
metric = 'varWithin'
dataType = 'abu'
log_tran='xy'
ylim=NULL
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------

results1 = empirVarBin
results2 = simVarBinFixed
results3 = simVarBinLogSer
metric = 'varWithin'
dataType = 'bin'
log_tran='xy'
ylim=NULL
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------

## compute residuals between observed and expected
resSorAbuFixed = getResid(empirSorAbu, simSorAbuFixed)
resSorAbuLogSer = getResid(empirSorAbu, simSorAbuLogSer)
resVarAbuFixed = getResid(empirVarAbu, simVarAbuFixed)
resVarAbuLogSer = getResid(empirVarAbu, simVarAbuLogSer)
##
resSorBinFixed = getResid(empirSorBin, simSorBinFixed)
resSorBinLogSer = getResid(empirSorBin, simSorBinLogSer)
resVarBinFixed = getResid(empirVarBin, simVarBinFixed)
resVarBinLogSer = getResid(empirVarBin, simVarBinLogSer)

plot(med.res ~ Dist, data=resSorAbuFixed, log='x')

dat = resSorAbuFixed
dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
sites = unique(dat$site)

avgresid = aggregate(avg.res ~ site, data=dat, function(x) sum(x^2) / length(x))
avgresid[order(avgresid[ , 2]), ]

shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
         "#FF8000")

par(mfrow=c(1,2))
plot(Metric.avg - avg.res ~ Metric.avg, data=dat,
     ylim=c(0, 1), type='n')
for(i in seq_along(sites)) {
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  tmp = subset(dat, site == sites[i])
  grains = unique(tmp$Comm)
  for(g in seq_along(grains)) {
    lines(Metric.avg - avg.res ~ Metric.avg, data=tmp, 
           subset = Comm == grains[g],
           pch=15, col=habcol[habindex], lwd=4)  
    lines(Metric.avg - exp.res ~ Metric.avg, data=tmp, 
           subset = Comm == grains[g], 
           pch=16, col=habcol[habindex], lwd=4)
  }  
}
abline(a=0, b=1)
##
h = hist(log2(dat$area), breaks=5, plot=F)
areaFac = cut(log2(dat$area), h$breaks)
areaFacUni = sort(unique(areaFac))
col = colorRampPalette(c('dodgerblue', 'red'))(length(areaFacUni))
plot(Metric.avg - avg.res ~ Metric.avg, data=dat,
     ylim=c(0, 1), type='n')
for(i in seq_along(sites)) {
  tmp = subset(dat, site == sites[i])
  grains = unique(tmp$area)
  for(g in seq_along(grains)) {
    areaFac = cut(log2(grains[g]), h$breaks)
    areaindex = match(areaFac, areaFacUni)
    lines(Metric.avg - avg.res ~ Metric.avg, data=tmp,
           subset = Comm == grains[g],
           pch=15, col=col[areaindex], lwd=4)  
    lines(Metric.avg - exp.res ~ Metric.avg, data=tmp,
           subset = Comm == grains[g],
           pch=16, col=col[areaindex], lwd=4)
  } 
}
abline(a=0, b=1)


par(mfrow=c(1,1))
plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim=c(-.5,.5))
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

## residuals vs distance
par(mfrow=c(1,1))
plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim=c(-.5,.5))
abline(h=0, lty=2, lwd=2)
for(i in seq_along(sites)) {
  tmp = subset(dat, site == sites[i])
  grains = unique(tmp$Comm)
  col = colorRampPalette(c('dodgerblue', 'red'))(length(grains))
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  lines(lowess(tmp$Dist, tmp$avg.res), col = habcol[habindex], lwd=2)
  lines(lowess(tmp$Dist, tmp$exp.res), col = habcol[habindex], lwd=2, lty=2)  
}

plot(avg.res ~ area, data=dat, log='x', ylim=c(-.5, .5))
for(i in seq_along(sites)) {
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  lines(avg.res ~ area, data=dat, subset= site == sites[i] ,
         col=habcol[habindex], pch=1, lwd=2)
  points(exp.res ~ area, data=dat, subset= site == sites[i],
         col=habcol[habindex], pch=1, lwd=2)  
}

## residuals vs area
plot(avg.res ~ area, data=dat, type='n', log='x', axes=F, frame.plot=F,
     xlim = c(0.1, 1e5), ylim=c(-.6, .6), xlab='', ylab='')
axis(side=1, cex.axis=1.75, padj=.5, lwd=8,
     at=10 ^ (-1:5))
axis(side=2, cex.axis=1.75, lwd=8)
abline(h=0, lwd=5)
for(i in seq_along(sites)) {
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  tmp = subset(dat, site == sites[i])
  lines(lowess(tmp$area, tmp$avg.res),
         col=habcol[habindex], pch=1, lwd=4)
  lines(lowess(tmp$area, tmp$exp.res), lty=2,
         col=habcol[habindex], pch=1, lwd=4)  
}

mk_legend('center', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')
mk_legend('center', hab, col=habcol, lty=3, lwd=7, cex=2, bty='n')


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
  
## compute summary statics for empricial and simualted results
empirStatsSorAbuMed = getStats(empirSorAbu, 'median')
empirStatsSorAbuAvg = getStats(empirSorAbu, 'average')

empirStatsSorAbuMed$ferp = empirStatsSorAbuMed$ferp[,,,-6]
empirStatsSorAbuAvg$ferp = empirStatsSorAbuAvg$ferp[,,,-6]
  
simStatsSorAbuFixedMed = getStats(simSorAbuFixed, 'median')
simStatsSorAbuFixedAvg = getStats(simSorAbuFixed, 'average')

simStatsSorAbuLogSerMed = getStats(simSorAbuLogSer, 'median')
simStatsSorAbuLogSerAvg = getStats(simSorAbuLogSer, 'average')

##
empirStatsVarAbuMed = getStats(empirVarAbu, 'median')
empirStatsVarAbuAvg = getStats(empirVarAbu, 'average')

empirStatsVarAbuMed$ferp = empirStatsVarAbuMed$ferp[,,,-6]
empirStatsVarAbuAvg$ferp = empirStatsVarAbuAvg$ferp[,,,-6]
  
simStatsVarAbuFixedMed = getStats(simVarAbuFixed, 'median')
simStatsVarAbuFixedAvg = getStats(simVarAbuFixed, 'average')

simStatsVarAbuLogSerMed = getStats(simVarAbuLogSer, 'median')
simStatsVarAbuLogSerAvg = getStats(simVarAbuLogSer, 'average')

##-----------------------------------------------------------------------------

empirStatsSorBinMed = getStats(empirSorBin, 'median')
empirStatsSorBinAvg = getStats(empirSorBin, 'average')

empirStatsSorBinMed$ferp = empirStatsSorBinMed$ferp[,,,-6]
empirStatsSorBinAvg$ferp = empirStatsSorBinAvg$ferp[,,,-6]
  
simStatsSorBinFixedMed = getStats(simSorBinFixed, 'median')
simStatsSorBinFixedAvg = getStats(simSorBinFixed, 'average')

simStatsSorBinLogSerMed = getStats(simSorBinLogSer, 'median')
simStatsSorBinLogSerAvg = getStats(simSorBinLogSer, 'average')

##
empirStatsVarBinMed = getStats(empirVarBin, 'median')
empirStatsVarBinAvg = getStats(empirVarBin, 'average')

empirStatsVarBinMed$ferp = empirStatsVarBinMed$ferp[,,,-6]
empirStatsVarBinAvg$ferp = empirStatsVarBinAvg$ferp[,,,-6]
  
simStatsVarBinFixedMed = getStats(simVarBinFixed, 'median')
simStatsVarBinFixedAvg = getStats(simVarBinFixed, 'average')

simStatsVarBinLogSerMed = getStats(simVarBinLogSer, 'median')
simStatsVarBinLogSerAvg = getStats(simVarBinLogSer, 'average')


## generate pdfs of r2 values, exponential or pwr function better
pl_r2 = function() {
  for( i in 1:3) {
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

stats = list(empirStatsSorAbuAvg, simStatsSorAbuFixedAvg, simStatsSorAbuLogSerAvg)
pdf('./figs/exp_vs_pwr_avg_sor_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsSorAbuMed, simStatsSorAbuFixedMed, simStatsSorAbuLogSerMed)
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

## plot summary statistics

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
  stats3 = simStatsSorAbuLogSerAvg
  stats1 = stats1[-grep('serp', names(stats1))]
  stats2 = stats2[-grep('serp', names(stats2))]
  stats3 = stats3[-grep('serp', names(stats3))]

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


  pl_coef(stats1, stats2, stats3)
  #
  cof = 'b1'
  pl_coef(stats1, stats2, stats3)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2, stats3)
  #
  cof = 'b1'
  pl_coef(stats1, stats2, stats3)
##------------------------------
  stats1 = empirStatsSorBinAvg
  stats2 = simStatsSorBinFixedAvg
  stats3 = simStatsSorBinLogSerAvg
  stats1 = stats1[-grep('serp', names(stats1))]
  stats2 = stats2[-grep('serp', names(stats2))]
  stats3 = stats3[-grep('serp', names(stats3))]

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
#####
  pl_coef(stats1, stats2, stats3)
  #
  cof = 'b1'
  pl_coef(stats1, stats2, stats3)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2, stats3)
  #
  cof = 'b1'
  pl_coef(stats1, stats2, stats3)
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

################ Abu Sorensen
pdf('empir_sim_sor_arith.pdf')
par(mfrow=c(4,4))
ylims = rbind(c(0.001,.4),c(.001,.05),c(.001,.05),c(.001,.4),
              c(.001,.7),c(.001,.05),c(.001,.05),c(.001,.05))
for (i in seq_along(shrtnames)) {
  results1 = eval(parse(text=paste('empirSorAbu$', shrtnames[i], sep='')))           
  results2 = eval(parse(text=paste('simSorAbuLogSer$', shrtnames[i], sep=''))) 
  results3 = eval(parse(text=paste('simSorAbuFixed$', shrtnames[i], sep=''))) 
  unigrains = unique(results1$Comm)
  plot(Metric.50 ~ Dist, data = results1, ylim=c(0, 1), type='n', main=shrtnames[i])
  for (j in seq_along(unigrains)) {
    lines(Metric.50 ~ Dist, data = results1, subset= Comm == unigrains[j], 
         col=j, lty=1, lwd=2, type='o')
    lines(Med ~ Dist, data = results3, subset= Comm == unigrains[j],
          col=j, lty=2, lwd=2, type='l')
    if (shrtnames[i] != 'cross')
      lines(Med ~ Dist, data = results2, subset= Comm == unigrains[j],
            col=j, lty=3, lwd=2, type='l')
  }
}

################ Binary Sorensen
par(mfrow=c(2,4))
results1 = empirSorBin
results2 = simSorBinLogSer
results3 = simSorBinFixed
ylims = rbind(c(0.001,.4),c(.001,.05),c(.001,.05),c(.001,.6),
              c(.3,1),c(.001,.05),c(.001,.05),c(.001,.05))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l',
       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}
dev.off()
################ Abu varWithin
par(mfrow=c(2,4))
results1 = empirVarAbu
results2 = simVarAbuLogSer
results3 = simVarAbuFixed
ylims = rbind(c(35,2000),c(.001,.05),c(.001,.05),c(.001,.4),
              c(.3,.7),c(.001,.05),c(.001,.05),c(.001,.05))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l')
#       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}

################ binary varWithin
par(mfrow=c(2,4))
results1 = empirVarBin
results2 = simVarBinLogSer
results3 = simVarBinFixed
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l')
#       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}

dev.off()
############################################ log -log
################ Abu Sorensen
pdf('empir_sim_sor_loglog.pdf')
par(mfrow=c(2,4))
results1 = empirSorAbu
results2 = simSorAbuLogSer
results3 = simSorAbuFixed
ylims = rbind(c(0.02,1),c(.001,1),c(.001,1),c(.02,.4),
              c(.1,1),c(.003,1),c(.003,1),c(.003,1))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='n',
       ylim=ylims[i,],main=shrtnames[i],log='xy')
  tmp = results2[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='palegreen')
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='green3',lwd=2,
        type='l')
  tmp = results3[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='dodgerblue')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='blue',
        lwd=2,type='l')
  lines(Exp ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  points(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,pch=19)
}

################ Binary Sorensen
par(mfrow=c(2,4))
results1 = empirSorBin
results2 = simSorBinLogSer
results3 = simSorBinFixed
ylims = rbind(c(0.05,1),c(.001,1),c(.001,1),c(.03,.8),
              c(.2,1),c(.003,1),c(.003,1),c(.003,1))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='n',
       ylim=ylims[i,],main=shrtnames[i],log='xy')
  tmp = results2[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='palegreen')
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='green3',lwd=2,
        type='l')
  tmp = results3[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='dodgerblue')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='blue',
        lwd=2,type='l')
  lines(Exp ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  points(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,pch=19)
}
dev.off()

par(mfrow=c(1,1))
plot(1:10,1:10,frame.plot=F,axes=F,xlab='',ylab='',type='n')
legend('center',c('Empirical','METE','METE fixed abu','RP'),col=c('black','green3',
       'blue','red'),cex=2,pch=c(19,rep(NA,3)),lty=c(NA,rep(1,3)),lwd=c(NA,rep(2,4)),bty='n')



################ Abu varWithin
par(mfrow=c(2,3))
results1 = empirVarAbuAvg
results2 = simVarAbuLogSer
results3 = simVarAbuFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(35,2000),main='bci',log='xy')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.1,1),main='cocoli',log='xy')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.1,4),main='sherman1',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(1,4),main='sherman2',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(1e3,6e4),main='serp',log='xy')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

################ binary varWithin
par(mfrow=c(2,3))
results1 = empirVarBinAvg
results2 = simVarBinLogSer
results3 = simVarBinFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(4,13),main='bci',log='xy')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.35,.55),main='cocoli',log='xy')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.6,1),main='sherman1',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(1,1.6),main='sherman2',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(1,4),main='serp',log='xy')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

dev.off()
