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

load('./sorensen/empirSorAbu_bisect.Rdata')
load('simulated_empirical_results_bisect.Rdata')

resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)
resSorAbuLogSer = get_ddr_resid(empirSorAbu, simSorAbuLogSer)


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

## altnerative plots for Fig. 2
## site specific plots

dat = resSorAbuFixed
dat = data.frame(dat, area = as.numeric(as.character(dat$Comm)))
sar_data = read.csv('./sar/empir_sars.csv')
sar_data$area = round(sar_data$area, 2)
dat = merge(dat, sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)
## subset so that has at least 20 individuals
dat = subset(dat, indiv >= 20)
## bring in habitat type
dat$hab = habitat[match(dat$site, shrtnm)]
## compute raw mete and RP values
dat$mete = dat$Metric.avg - dat$avg.res
dat$rp = dat$Metric.avg - dat$exp.res

sites = unique(dat$site)

site_names = as.character(as.matrix(read.table('./data/shrtnames.txt'))) 
site_names = "bci, sherman1, cocoli1, luquillo, bryan, bigoak, oosting, rocky, bormann, woodbridge, baldmnt, landsend, graveyard, ferp, serp, cross"
site_names = unlist(strsplit(site_names, split=', '))
site_titles = sub('1', '', site_names)

capwords = function(s, strict = FALSE) {
  cap = function(s) paste(toupper(substring(s, 1, 1)),
{s = substring(s, 2); if(strict) tolower(s) else s},
                          sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

site_titles = capwords(site_titles)
site_titles[1] = "BCI"
site_titles[11] = "Bald Mtn."
site_titles[14] = "UCSC"
site_titles[15] = "Serpentine"
site_titles[16] = "Cross Timbers"

col = c('red', 'dodgerblue')
lty = c(2, 3)

## arith-arith plots
png('./figs/ddr_arith_resid_by_sites.png', width=480 * 4, height=480 * 4)
  par(mfrow=c(4,4))
  for (i in seq_along(site_names)) {
    metrics = c('avg.res', 'exp.res')
    true = as.character(dat$site) == site_names[i]
    xlim = (range(dat$Dist[true], na.rm=T))
    ylim = (range(dat[true , metrics, ], na.rm=T))
    xlim = c(floor(xlim[1]), ceiling(xlim[2]))
    plot((Dist) ~ (Metric.avg), data=dat[true, ],
         xlim=xlim, ylim=ylim, type='n',
         xlab='', ylab='', frame.plot=F, axes=F)
    addAxis(side=1, cex.axis=3, padj=.75)
    addAxis(side=2, cex.axis=3, padj=0)
    mtext(side=3, paste(site_titles[i], '-', unique(dat$hab[true])),
          cex=2)
    mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
    abline(h=0, col='grey', lwd=5)
    for (j in seq_along(metrics)) {
      lines(lowess((dat[true, 'Dist']), 
                   (dat[true, metrics[j]]), f=3/4), lwd=5, lty=lty[j], col=col[j])
    }
    if(i == 1)
      legend('right', 
             c('observed', 'recursive, METE SAD', 'recursive, observed SAD',
               'non-recursive, METE SAD', 'non-recursive, observed SAD'),
             pch=c(1, rep(NA, 4)), col=c(1, col), cex=2.5, bty='n',
             lwd=c(3, rep(5, 4)), lty=c(NA, lty))
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