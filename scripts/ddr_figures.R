source('./scripts/spat_functions.R')

## Figure 2: 4 x 4 panel site specific DDR graphic------------------------------

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

fixed = resSorAbuFixed
logser = resSorAbuLogSer
fixed = data.frame(fixed, area = as.numeric(as.character(fixed$Comm)))
logser = data.frame(logser, area = as.numeric(as.character(logser$Comm)))
sar_data = read.csv('./sar/empir_sars.csv')
sar_data$area = round(sar_data$area, 2)
fixed = merge(fixed, sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)
logser = merge(logser,  sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)

## bring in habitat type
fixed$hab = habitat[match(fixed$site, shrtnm)]
logser$hab =habitat[match(logser$site, shrtnm)] 

## compute raw METE and RPM values
fixed$mete = fixed$Metric.avg - fixed$avg.res
fixed$rp = fixed$Metric.avg - fixed$exp.res
logser$mete = logser$Metric.avg - logser$avg.res
logser$rp = logser$Metric.avg - logser$exp.res

sites = unique(fixed$site)

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

col = c('red', 'dodgerblue', 'black')
lty = c(1, 3)

find_middle = function(len) floor(len/2) + (len %%2)

shrtnames = read.table('./data/shrtnames.txt')
shrtnames = as.character(as.matrix(shrtnames[1 , ]))
bisect = read.table('./data/bisect.txt')
grain_fine = read.table('./data/grain_fine.txt')
A0 = 2^bisect[ , 2] * grain_fine
bisect_num = ifelse(bisect[ , 2] %% 2 == 0, 8, 9)
area_of_interest = as.numeric(round(A0 / 2 ^ bisect_num, 2))

png('./figs/ddr_arith_by_sites.png', width=480 * 4, height=480 * 4)
par(mfrow=c(4,4))
for (i in seq_along(site_names)) {
  metrics = c('mete', 'rp', 'Metric.avg')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  index = match(site_names[i], shrtnames)
  true_fix = true_fix & fixed$area == area_of_interest[index]
  true_log = true_log & logser$area == area_of_interest[index]
  xlim = (range(fixed$Dist[true_fix], na.rm=T))
  if (sum(true_log) > 0)
    ylim = (range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = (range(fixed[true_fix , metrics], na.rm=T))
  x = ifelse(ylim[1] < 0, -1, 1)
  ylim = c(floor((ylim[1] %% x) * 10) / 10, 
           ceiling((ylim[2] %% 1) * 10) / 10)
  plot(Metric.avg ~ Dist, data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='o', lwd=3, cex=2, pch=19,
       xlab='', ylab='', frame.plot=F, axes=F, log='')
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0)
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  metrics = metrics[-3]
  for (j in seq_along(metrics)) {
    lines(fixed[true_fix, 'Dist'], fixed[true_fix, metrics[j]],
          lwd=3, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 1)
    legend('right', 
           c('observed', 'recursive, observed SAD',
             'recursive, METE SAD','random, observed SAD'),
           col=c('black', rep(col, each=2)), cex=2.5, bty='n',
           lwd=rep(5, 5), lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()

png('./figs/ddr_loglog_by_sites.png', width=480 * 4, height=480 * 4)
par(mfrow=c(4,4))
for (i in seq_along(site_names)) {
  metrics = c('mete', 'rp', 'Metric.avg')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  index = match(site_names[i], shrtnames)
  true_fix = true_fix & fixed$area == area_of_interest[index]
  true_log = true_log & logser$area == area_of_interest[index]
  xlim = range(fixed$Dist[true_fix], na.rm=T)
  xliml2= log2(xlim)
  xends = c(floor(xliml2[1]), ceiling(xliml2[2]))
  xlim = 2^xends
  if (sum(true_log) > 0)
    ylim = (range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = (range(fixed[true_fix , metrics], na.rm=T))
  yliml2= log2(ylim)
  yends = c(floor(yliml2[1]), ceiling(yliml2[2]))
  ylim = 2^yends
  plot(Metric.avg ~ Dist, data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='o', lwd=3, cex=2, pch=19,
       xlab='', ylab='', frame.plot=F, axes=F, log='xy')
  xticks = 2^(xends[1]:xends[2]) 
  yticks = 2^(yends[1]:yends[2])
  addAxis(side=1, cex.axis=3, padj=.75, at=xticks, lab=as.character(xticks))
  addAxis(side=2, cex.axis=3, padj=0, at=yticks, lab=as.character(yticks))
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  metrics = metrics[-3]
  for (j in seq_along(metrics)) {
    lines(fixed[true_fix, 'Dist'], fixed[true_fix, metrics[j]],
          lwd=3, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 1)
    legend('bottomleft', 
            c('observed', 'recursive, observed SAD',
              'recursive, METE SAD','random, observed SAD'),
            col=c('black', rep(col, each=2)), cex=2.5, bty='n',
            lwd=rep(5, 5), lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()


png('./figs/ddr_log2log2_by_sites.png', width=480 * 4, height=480 * 4)
par(mfrow=c(4,4))
for (i in seq_along(site_names)) {
  metrics = c('mete', 'rp', 'Metric.avg')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  index = match(site_names[i], shrtnames)
  true_fix = true_fix & fixed$area == area_of_interest[index]
  true_log = true_log & logser$area == area_of_interest[index]
  xlim = log2(range(fixed$Dist[true_fix], na.rm=T))
  if (sum(true_log) > 0)
    ylim = log2(range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = log2(range(fixed[true_fix , metrics], na.rm=T))
  xlim = c(floor(xlim[1]), ceiling(xlim[2]))
  ylim = c(floor(ylim[1]), ceiling(ylim[2]))
  plot(log2(Metric.avg) ~ log2(Dist), data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='o', lwd=3, cex=2, pch=19,
       xlab='', ylab='', frame.plot=F, axes=F)
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0, at=ylim[1]:ylim[2])
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  metrics = metrics[-3]
  for (j in seq_along(metrics)) {
    lines(log2(fixed[true_fix, 'Dist']), log2(fixed[true_fix, metrics[j]]),
          lwd=3, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(log2(logser[true_log, 'Dist']), log2(logser[true_log, metrics[j]]),
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 1)
    legend('bottomleft', 
           c('observed', 'recursive, observed SAD',
             'recursive, METE SAD','random, observed SAD'),
           col=c('black', rep(col, each=2)), cex=2.5, bty='n',
           lwd=rep(5, 5), lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()

## Figure 3: one-to-one plots---------------------------------------------------
titles = c('recursive, METE SAD', 'recursive, observed SAD',
           'random, observed SAD')

png('./figs/ddr_one_to_one.png',
    width = 480 * 3, height= 480 * 1)
par(mfrow=c(1,3))
lims = c(0, .85)
for(i in 1:3) {
  metrics = c('mete', 'mete', 'rp')
  if (i == 1) {
    x = logser[ , metrics[i]]
    y = logser[ , 'Metric.avg']
  }
  else {
    x = fixed[ , metrics[i]]
    y = fixed[ , 'Metric.avg']
  }
  plot((x), (y), type='n', axes=F, frame.plot=F, xlab='', ylab='',
       xlim=lims, ylim=lims)
  if (i == 1) {
    true = logser$indiv >= 20
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  else {
    true = fixed$indiv >= 20
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  addAxis(side=1)
  addAxis(side=2)
  lines(lims, lims, lwd=2)
  mtext(side=3, titles[i], cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
}
dev.off()

png('./figs/ddr_loglog_one_to_one.png',
    width = 480 * 3, height= 480 * 1)
par(mfrow=c(1,3))
lims = c(2^-7, 1)
for(i in 1:3) {
  metrics = c('mete', 'mete', 'rp')
  if (i == 1) {
    x = logser[ , metrics[i]]
    y = logser[ , 'Metric.avg']
  }
  else {
    x = fixed[ , metrics[i]]
    y = fixed[ , 'Metric.avg']
  }
  plot((x), (y), type='n', axes=F, frame.plot=F, xlab='', ylab='',
       xlim=lims, ylim=lims, log='xy')
  if (i == 1) {
    true = logser$indiv >= 20
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  else {
    true = fixed$indiv >= 20
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  ticks = round(2^(-7:0), 3)
  addAxis(side=1, at = ticks, lab = as.character(ticks))
  addAxis(side=2, at = ticks, lab = as.character(ticks), padj=0)
  lines(lims, lims, lwd=2)
  mtext(side=3, titles[i], cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
}
dev.off()

png('./figs/ddr_log2log2_one_to_one.png',
    width = 480 * 3, height= 480 * 1)
  par(mfrow=c(1,3))
  lims = c(-7, 0)
  for(i in 1:3) {
    metrics = c('mete', 'mete', 'rp')
    if (i == 1) {
      x = logser[ , metrics[i]]
      y = logser[ , 'Metric.avg']
    }
    else {
      x = fixed[ , metrics[i]]
      y = fixed[ , 'Metric.avg']
    }
    plot(log2(x), log2(y), type='n', axes=F, frame.plot=F, xlab='', ylab='',
         xlim=lims, ylim=lims)
    if (i == 1) {
      true = logser$indiv >= 20
      points(log2(x[true]), log2(y[true]), pch=19)
      points(log2(x[!true]), log2(y[!true]), pch=19, col='grey')
    } 
    else {
      true = fixed$indiv >= 20
      points(log2(x[true]), log2(y[true]), pch=19)
      points(log2(x[!true]), log2(y[!true]), pch=19, col='grey')
    }
    ticks = lims[1]:lims[2]
    addAxis(side=1, at = ticks, lab = as.character(ticks))
    addAxis(side=2, at = ticks, lab = as.character(ticks), padj=0)
    lines(lims, lims, lwd=2)
    mtext(side=3, titles[i], cex=2)
    mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  }
dev.off()


## Supplemental site speciefic arith-arith residual plots -----------------------------------
png('./figs/ddr_arith_resid_by_sites.png', width=480 * 4, height=480 * 4)
par(mfrow=c(4,4))
for (i in seq_along(site_names)) {
  metrics = c('avg.res', 'exp.res')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  area = sort(unique(fixed$area[true_fix]))
  area = area[find_middle(length(area))]
  true_fix = true_fix & fixed$area == area
  true_log = true_log & logser$area == area
  xlim = (range(fixed$Dist[true_fix], na.rm=T))
  if (sum(true_log) > 0)
    ylim = (range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = (range(fixed[true_fix , metrics], na.rm=T))
  xlim = c(floor(xlim[1]), ceiling(xlim[2]))
  #ylim = round(ylim, 1)
  x = ifelse(ylim[1] < 0, -1, 1)
  ylim = c(floor((ylim[1] %% x) * 10) / 10, 
           ceiling((ylim[2] %% 1) * 10) / 10)
  plot((Dist) ~ (Metric.avg), data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='n',
       xlab='', ylab='', frame.plot=F, axes=F)
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0)
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  abline(h=0, col='grey', lwd=5)
  for (j in seq_along(metrics)) {
    lines(fixed[true_fix, 'Dist'], fixed[true_fix, metrics[j]],
          lwd=3, lty=lty[1], col=col[j], type='o')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='o')
    }  
  }
  if(i == 1)
    legend('right', 
           c('recursive, observed SAD', 'recursive, METE SAD',
             'random, observed SAD'),
           col=rep(col, each=2), cex=2.5, bty='n', lwd=rep(5, 4), lty=lty)
}
dev.off()


## alternative loess approach with SE


## log-arith lowess plots
tiff('./figs/ddr_log_resid_by_sites.tiff', width=480 * 4, height=480 * 4)
par(mfrow=c(4,4))
for (i in seq_along(site_names)) {
  metrics = c('avg.res', 'exp.res')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  xlim = (range(fixed$Dist[true_fix], na.rm=T))
  if (sum(true_log) > 0)
    ylim = (range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = (range(fixed[true_fix , metrics], na.rm=T))
  xlim = c(floor(xlim[1]), ceiling(xlim[2]))
  xlim[1] = ifelse(xlim[1] == 0, 1, xlim[1])
  #ylim = round(ylim, 1)
  x = ifelse(ylim[1] < 0, -1, 1)
  ylim = c(floor((ylim[1] %% x) * 10) / 10, 
           ceiling((ylim[2] %% 1) * 10) / 10)
  plot((Dist) ~ (Metric.avg), data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='n',
       xlab='', ylab='', frame.plot=F, axes=F, log='x')
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0)
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  abline(h=0, col='grey', lwd=5)
  for (j in seq_along(metrics)) {
    lines(lowess((fixed[true_fix, 'Dist']), 
                 (fixed[true_fix, metrics[j]]), f=7/8), lwd=3, lty=lty[1], col=col[j])
    if (j == 1 & sum(true_log) > 0) ## logser results only for METE model
      lines(lowess((logser[true_log, 'Dist']), 
                   (logser[true_log, metrics[j]]), f=7/8), lwd=5, lty=lty[2], col=col[j])
    
  }
  if(i == 1)
    legend('right', 
           c('recursive, observed SAD', 'recursive, METE SAD',
             'random, observed SAD'),
           col=rep(col, each=2), cex=2.5, bty='n', lwd=rep(5, 4), lty=lty)
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