setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

##------------------------------------------------------------------------------
## Figure 2 - example Simulated DDR with simulation results
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
          quants=F, col=col, lwd=5, add=TRUE, type='p', cex=1.5, pch=19)
for (g in seq_along(grains)) {
  mod = lm(log(Avg) ~ log(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  lines(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]], exp(predict(mod)),
        col=col[g], lwd=linelwd)  
}

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
###

col = colorRampPalette(c('dodgerblue', 'red'))(5)
ddr$ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg) / grains
xlab = 'log(N/S) / log(Nbar/Sbar)'
plot(10^b0 ~ ratio , data=ddr, xlab='', ylab='', type='n',
     frame.plot=F, axes=F, ylim=c(0, 1), xlim=c(0,80))
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
     frame.plot=F, axes=F, xlim=c(0, 80))
axis(side=1, cex.axis=1.75, padj=.5, lwd=axislwd)
axis(side=2, cex.axis=1.75, lwd=axislwd) 
for (g in seq_along(grains)) {
  tmp = subset(ddr, grains == grains[g])
  lines(lowess(tmp$ratio, tmp$b1, f= 1/3), col=col[g], lwd=linelwd)
}  

mk_legend('center', legend=grains, col=col,
          lty=1, bty='n', cex = 2, lwd=8)

##----------------------------------------------------------------------------
## empirical DDR pattern at a single scale

## start with a good fit for slope by METE
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

##----------------------------------------------------------------------------
## Supplemental Figure 1
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

##----------------------------------------------------------------------------
## emprical DDR residuals
load('./sorensen/empirSorAbu.Rdata') 
## drop ferp last grain
empirSorAbu$ferp = empirSorAbu$ferp[-nrow(empirSorAbu$ferp),]
resSorAbuFixed = getResid(empirSorAbu, simSorAbuFixed)
resSorAbuLogSer = getResid(empirSorAbu, simSorAbuLogSer)

dat = resSorAbuFixed
dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
sites = unique(dat$site)

shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
         "#FF8000")

## residuals vs distance
par(mfrow=c(1,1))
plot(avg.res ~ Dist, data=dat, log='x', axes=F, frame.plot=F,
     type='n', ylim=c(-.5,.5), xlab='', ylab='')
axis(side=1, cex.axis=1.75, padj=.5, lwd=8,
     at= .5 * 10^(0:3))
axis(side=2, cex.axis=1.75, lwd=8)
abline(h=0, lwd=6)
for(i in seq_along(sites)) {
  tmp = subset(dat, site == sites[i])
  grains = unique(tmp$Comm)
  col = colorRampPalette(c('dodgerblue', 'red'))(length(grains))
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  lines(lowess(tmp$Dist, tmp$avg.res), col = habcol[habindex], lwd=4)
  lines(lowess(tmp$Dist, tmp$exp.res), col = habcol[habindex], lwd=4, lty=2)  
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

##----------------------------------------------------------------------------empirStatsSorAbuMed = getStats(empirSorAbu, 'median')
## METE DDR scale colapse with data
empirStatsSorAbuAvg = getStats(empirSorAbu, 'average')
empirStatsSorAbuAvg$ferp = empirStatsSorAbuAvg$ferp[,,,-6]

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

##----------------------------------------------------------------------------
## empirical SAR pattern

fileNames = dir('./sar')

empirFiles = grep('empir_sar.csv', fileNames)
empir = vector('list', length(empirFiles))
names(empir) = sub('_empir_sar.csv', '', fileNames[empirFiles])
for (i in seq_along(empirFiles))
  empir[[i]] = read.csv(paste('./sar/', fileNames[empirFiles[i]], sep=''))

meteFiles = grep('mete_sar.txt', fileNames)
mete = vector('list', length(meteFiles))
names(mete) = sub('_mete_sar.txt', '', fileNames[meteFiles])
for (i in seq_along(meteFiles)) {
  mete[[i]] = read.csv(paste('./sar/', fileNames[meteFiles[i]], sep=''))
  mete[[i]]$area = mete[[i]]$area * empir[[i]]$area[1]
}

meteavgFiles = grep('mete_sar_avgs.csv', fileNames)
meteavg = vector('list', length(meteavgFiles))
names(meteavg) = sub('_mete_sar_avgs.csv', '', fileNames[meteavgFiles])
for (i in seq_along(meteavgFiles)) {
  meteavg[[i]] = read.csv(paste('./sar/', fileNames[meteavgFiles[i]], sep=''))
  meteavg[[i]] = meteavg[[i]][ , -1]
  Amin = empir[[match(names(meteavg)[i], names(empir))]]$area[1]
  meteavg[[i]]$grains = meteavg[[i]]$grains * Amin
}

## average cocoli and sherman plots
empir$sherman1 = (empir$sherman1 + empir$sherman2) / 2
empir$cocoli1 = (empir$cocoli1 + empir$cocoli2) / 2
mete$sherman1 = (mete$sherman1 + mete$sherman2) / 2
mete$cocoli1 = (mete$cocoli1 + mete$cocoli2) / 2
meteavg$sherman1 = (meteavg$sherman1 + meteavg$sherman2) / 2
meteavg$cocoli1 = (meteavg$cocoli1 + meteavg$cocoli2) / 2
empir = empir[-match(c('sherman2','sherman3','cocoli2'), names(empir))]
mete = mete[-match(c('sherman2','sherman3','cocoli2'), names(mete))]
meteavg = meteavg[-match(c('sherman2','cocoli2'), names(meteavg))]

load('./sar/expected_empir_sars.Rdata')

addCI = function(x, y.lo, y.hi, col, data=NULL) {
  if (!is.null(data)) {
    x = eval(parse(text=paste(data, '$', x, sep='')))
    y.lo = eval(parse(text=paste(data, '$', y.lo, sep='')))
    y.hi = eval(parse(text=paste(data, '$', y.hi, sep='')))
  }  
  xvals = c(x, rev(x))
  yvals = c(y.lo, rev(y.hi))
  polygon(xvals, yvals, border=NA, col=col)
}

site = 'bigoak'
purple = rgb(112, 48, 160, maxColorValue=160)
lightblue = "#1AB2FF"

i = match(site, names(mete))
    plot(sr ~ area, data=mete[[i]], log='xy',
         ylim=c(0.2, 52),
         xlim=c(1, 1.5e4),
         type='n', xlab='', ylab='',frame.plot=F, axes=F)
    axis(side=1, cex.axis=1.75, padj=.5, lwd=8)
    ticks = 0.2 * 2^seq(0, 8, 2)
#    ticks[ticks > 1] = round(ticks[ticks >1])
    axis(side=2, cex.axis=1.75, lwd=8,
         at = ticks) 
    ## mete CI
    dat = meteavg[[match(names(mete)[i], names(meteavg))]]
    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
    lines(sr.avg ~ grains, data=dat, col=lightblue, lwd=4, lty=2)
    ## RP CI
    dat = srExp[[match(names(mete)[i], names(srExp))]]
#    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='pink')
    lines(srCol~ grains, data=dat, col=purple, lwd=4, lty=2)
    ## analytical mete    
#    lines(sr ~ area, data=mete[[i]], type='o')
    ## data
    lines(richness ~ area, data = empir[[i]], pch=19, type='o',
          lwd=4, cex=1.25)

##----------------------------------------------------------------------------
## empirical SAR residuals

for (i in seq_along(empir)) {
  site = names(empir)[i]
  if (site == 'cross') {
    index = match(site, names(mete))
    mete.sr = mete[[index]]$sr
    area = mete[[index]]$area
    obs.sr = empir[[i]]$richness[match(area, empir[[i]]$area)]
    mete.res = obs.sr - mete.sr
    index = match(site, names(srExp))
    rp.sr = srExp[[index]]$srCol[match(area, srExp[[index]]$grains)]
    rp.res = obs.sr - rp.sr
  }
  else {
    area = empir[[i]]$area
    obs.sr = empir[[i]]$richness
    index = match(site, names(meteavg))
    mete.sr =  meteavg[[index]]$sr.avg
    mete.res = obs.sr - mete.sr
    index = match(site, names(srExp))
    rp.sr = srExp[[index]]$srCol
    rp.res = obs.sr - rp.sr
  }
  if (exists('sr_res'))
    sr_res = rbind(sr_res,
                   data.frame(site, area, obs.sr, mete.sr, rp.sr, mete.res, rp.res))
  else
    sr_res = data.frame(site, area, obs.sr, mete.sr, rp.sr, mete.res, rp.res)
}

sites = unique(sr_res$site)
plot(mete.res ~ area, data=sr_res, log='x', ylim=c(-40, 40), 
     type='n', frame.plot=F, axes=F, xlab='', ylab='',
     xlim = c(0.1, 1e6))
axis(side=1, cex.axis=1.75, padj=.5, lwd=8,
     at=10 ^ (-1:6))
axis(side=2, cex.axis=1.75, lwd=8)
abline(h=0, lwd=5)
for (i in seq_along(sites)) {
  habindex = match(habitat[match(sites[i], shrtnm)], hab)
  lines(mete.res ~ area, data=sr_res, subset= site == sites[i],
        lwd=4, col=habcol[habindex])
  lines(rp.res ~ area, data=sr_res, subset= site == sites[i],
        lwd=4, col=habcol[habindex], lty=2)
}





