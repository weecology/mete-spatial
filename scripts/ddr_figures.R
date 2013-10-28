setwd('~/maxent/spat')
source('./scripts/spat_functions.R')

## Figure 2: 3 panel DDR graphic-------------------------------------------------
## A) example DDR, B) DDR METE residuals, C) DDR RP residuals
##
## set up graphic parameters 
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak woodland', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
        "#FF8000")

load('./sorensen/empirSorAbu.Rdata')
load('simulated_empirical_results.Rdata')
resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)

#pdf('./figs/fig2_ddr_empirical.pdf', width= 7 * 2, height= 7)
#windows(width= 7 * 3, height= 7 * 2)

  par(mfrow=c(1, 3))

  ## panel A: example empirical DDR with predictions
  tmp = empirSorAbu['bigoak']
  obs = empirSorAbu$'bigoak'
  exp = simSorAbuLogSer$'bigoak'
  grains = unique(obs$Comm)
  plot(Metric.avg ~ Dist, data=obs,
       ylim=c(0.02, .64), xlim=c(2,128), type='n', frame.plot=F, axes=F,
       xlab='', ylab='', log='xy')
  addAxis(1, at = 2^(1:7))
  addAxis(2, at= 0.02 * 2^(0:5), padj=-.5)
  g = 2
  ## add mete  
  true = exp$Comm == grains[g]
  tmpexp = exp[true,]
  lines(Avg ~ Dist, data=tmpexp, col=1, lty=1, lwd=3)
  ## RP
  lines(Exp.avg ~ Dist, data=obs, subset=Comm == grains[g], lty=1,
        col='grey', lwd=3)  
  ## data
  lines(Metric.avg ~ Dist, data=obs, subset=Comm == grains[g],
        lty=1, lwd=3, col=1, type='p', pch=19, cex=2)
  legend('topright', c('Observed', 'METE', 'RP'), pch=c(19, NA, NA), lwd=c(NA,5,5),
         col=c(1, 1, 'grey'), cex=2, bty='n')  
  ## panels B & C: emprical DDR residuals
  dat = resSorAbuFixed
  dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
  sar_data = read.csv('./sar/empir_sars.csv')
  sar_data$area = round(sar_data$area, 2)
  dat = merge(dat, sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)
  ## subset so that has at least 20 individuals
  dat = subset(dat, indiv >= 20)
  ## bring in habitat type
  dat$hab = habitat[match(dat$site, shrtnm)]
  sites = unique(dat$site)
  for(j in 1:2){
    if(j == 1)
      main = 'METE'
    else
      main = 'RP'
    plot(avg.res / Metric.avg ~ Dist, data=dat, log='x', type='n', ylim = c(-.1, .8), xlim=c(.5,512),
         xlab='', ylab='', axes=F, frame.plot=F)
    addAxis(1, at=2^seq(-1, 9, 2))
    addAxis(2, padj=-.5)
    abline(h=0, lwd=4, lty=2)
    for(i in seq_along(sites)) {
      tmp = subset(dat, site == sites[i])
      grains = unique(tmp$Comm)
      habindex = match(habitat[match(sites[i], shrtnm)], hab)
      if(j == 1) {
        lines(lowess(tmp$Dist, tmp$avg.res / tmp$Metric.avg, f=.75), col = habcol[habindex], lwd=3)
      }
      else {
        lines(lowess(tmp$Dist, tmp$exp.res / tmp$Metric.avg, f=.75), col = habcol[habindex], lwd=3)   
      }        
    }  
  }  
dev.off()

## Supplemental Figure 1--------------------------------------------------------
## r2 of model fits to simulated results
load('./sorensen/simSorAbuAvg.Rdata')
S = round(10^seq(log10(10), log10(100), length.out=20))
N = round(10^seq(log10(120), log10(5e5), length.out=20))
#stats = getSimStats(simSorAbuAvg, S, N)

#pdf('./figs/sup_fig1_r2_sim_&_empir_ddr.pdf', width = 7 * 2, height= 7)
#windows(width= 7 * 2, height=7)

meth='wtr'
dpwr = density(stats['pwr', 'r2', meth, , , ], na.rm = TRUE)
hpwr = hist(stats['pwr', 'r2', meth, , , ], plot=F)
dexp = density(stats['exp', 'r2', meth, , , ], na.rm = TRUE)
hexp = hist(stats['exp', 'r2', meth, , , ], plot=F)

xexp = dexp$x
yexp = dexp$y / sum(hexp$density) * dexp$n

xpwr = c(min(xexp), dpwr$x)
ypwr = c(0, dpwr$y / sum(hpwr$density) * dpwr$n)

xlims = range(c(xpwr, xexp, 1))
ylims = range(c(ypwr, yexp))

linelwd = 3

par(mfrow=c(1,2))

plot(xpwr, ypwr, type='l', lty=3, lwd=linelwd, xlim=round(xlims,1), ylim=ylims, col='black',
     xlab='', ylab='', frame.plot=F, axes=F)
mtext(side=3, 'METE parameter space', cex=2)
addAxis1(at=c(.7, .8, .9, 1))
addAxis2()
addxlab(expression('Coefficient of Determination, ' * italic(R^2)),
        padj=2)
addylab('Freqency')
lines(xpwr, ypwr,  lwd=linewd, col='black')
lines(xexp, yexp, lwd=linelwd, col='grey')

##
load('./sorensen/empirSorAbu.Rdata') 
empir_stats = getStats(empirSorAbu, 'average')

r2pwr = unlist(sapply(empir_stats, function(x) x['pwr', 'r2', 'wtr',]))
r2exp = unlist(sapply(empir_stats, function(x) x['exp', 'r2', 'wtr',]))
dpwr = density(r2pwr, na.rm=T)
hpwr = hist(r2pwr, plot=F)
dexp = density(r2exp, na.rm=T)
hexp = hist(r2exp, plot=F)

xexp = dexp$x
yexp = dexp$y / sum(hexp$density) * dexp$n

xpwr = c(min(xexp), dpwr$x)
ypwr = c(0, dpwr$y / sum(hpwr$density) * dpwr$n)


plot(xpwr, ypwr, type='n', xlim=range(c(xpwr, xexp)), ylim=range(c(ypwr, yexp)),
     xlab='', ylab='', frame.plot=F, axes=F)
mtext(side=3, 'Empirical datasets', cex=2)
addAxes()
addxlab(expression('Coefficient of Determination, ' * italic(R^2)),
        padj=2)
addylab('Freqency')
lines(xpwr, ypwr,  lwd=linelwd, col='black')
lines(xexp, yexp, lwd=linelwd, col='grey')

legend('topleft', c('Exponential Model', 'Power  Model'),
       lty=1, lwd=6, bty='n', cex=1.5, col=c( 'grey', 'black'))


dev.off()

## Supplemental Figure 2--------------------------------------------------------

source('./scripts/spat_sar_load_and_avg_data.R')

## set up graphic parameters 

shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak woodland', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
           "#FF8000")

resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)
dat = resSorAbuFixed
dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
sar_data = read.csv('./sar/empir_sars.csv')
sar_data$area = round(sar_data$area, 2)
dat = merge(dat, sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)
## subset so that has at least 20 individuals
dat = subset(dat, indiv >= 20)
## bring in habitat type
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
dat$hab = habitat[match(dat$site, shrtnm)]
sites = unique(dat$site)
#pdf('./figs/sup_fig2_raw_ddr_residuals.pdf', width=7 * 2, height=7)
#windows(width=7 * 3, height=7 * 2)
## raw emprical DDR residuals
par(mfrow=c(1,2))
for(j in 1:2){
  if(j == 1)
    main = 'METE'
  else
  main = 'RP'
  plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim = c(-.1, .6), xlim=c(.5,512),
       xlab='', ylab='', axes=F, frame.plot=F)
  addAxis(1, at=2^seq(-1, 9, 2))
  addAxis(2)
  abline(h=0, lwd=4, lty=2)
  addxlab('Distance (m)', padj=2.25)
  mtext(side=2,'Raw Residuals',cex=2,padj=-1.5)
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
dev.off()