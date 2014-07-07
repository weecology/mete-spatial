print('Generating figures, ...')

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
fixed$mete.50 = fixed$Metric.50 - fixed$med.res
fixed$rp = fixed$Metric.avg - fixed$exp.res
logser$mete = logser$Metric.avg - logser$avg.res
logser$mete.50 = logser$Metric.50 - logser$med.res
logser$rp = logser$Metric.avg - logser$exp.res

site_names = as.character(as.matrix(read.table('./data/shrtnames.txt'))) 
site_titles = site_names

if (length(site_names) >= 16) {
  site_names = "bci, sherman1, cocoli1, luquillo, bryan, bigoak, oosting, rocky, bormann, woodbridge, baldmnt, landsend, graveyard, ferp, serp, cross"
  site_names = unlist(strsplit(site_names, split=', '))
  site_titles = sub('1', '', site_names)
  site_titles = capwords(site_titles)
  site_titles[1] = "BCI"
  site_titles[11] = "Bald Mtn."
  site_titles[14] = "UCSC"
  site_titles[15] = "Serpentine"
  site_titles[16] = "Cross Timbers"
}

col = c('red3', 'lightskyblue3', 'black')
lty = c(1, 3)

find_middle = function(len) floor(len/2) + (len %%2)

shrtnames = read.table('./data/shrtnames.txt')
shrtnames = as.character(as.matrix(shrtnames[1 , ]))
bisect_fine = read.table('./data/bisect_fine.txt')
grain_fine = read.table('./data/grain_fine.txt')
A0 = 2^bisect_fine * grain_fine
bisect_num = ifelse(bisect_fine %% 2 == 0, 8, 9)
area_of_interest = as.numeric(round(A0 / 2 ^ bisect_num, 2))

if (length(site_names) == 2) {
  pltpar = list()
  pltpar$wd_mult = 2
  pltpar$ht_mult = 1
  pltpar$mfrow = c(1, 2)
}
if (length(site_names) == 16) {
  pltpar = list()
  pltpar$wd_mult = 4.25
  pltpar$ht_mult = 4
  pltpar$mfrow = c(4, 4)
}

png('./figs/ddr_arith_by_sites.png',
    width=480*pltpar$wd_mult, height=480*pltpar$ht_mult)
par(mfrow=pltpar$mfrow)
for (i in seq_along(site_names)) {
  metrics = c('mete', 'rp', 'Metric.avg')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  index = match(site_names[i], shrtnames)
  true_fix = true_fix & fixed$area == area_of_interest[index]
  true_log = true_log & logser$area == area_of_interest[index]
  xlim = (range(c(0,fixed$Dist[true_fix]), na.rm=T))
  if (sum(true_log) > 0)
    ylim = (range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = (range(fixed[true_fix , metrics], na.rm=T))
  x = ifelse(ylim[1] < 0, -1, 1)
  ylim = c(floor((ylim[1] %% x) * 10) / 10, 
           ceiling((ylim[2] %% 1) * 10) / 10)
  plot(Metric.avg ~ Dist, data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='o', lwd=5, cex=2, pch=19,
       xlab='', ylab='', frame.plot=F, axes=F, log='')
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0)
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  metrics = metrics[-3]
  for (j in seq_along(metrics)) {
    lines(fixed[true_fix, 'Dist'], fixed[true_fix, metrics[j]],
          lwd=5, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 13)
    legend('bottomleft', 
           c('observed', 'recursive, observed SAD',
             'recursive, METE SAD','random, observed SAD'),
           col=c('black', rep(col, each=2)), cex=3, bty='n',
           lwd=8, lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()

png('./figs/ddr_arith_by_sites_median.png',
    width=480*pltpar$wd_mult, height=480*pltpar$ht_mult)
par(mfrow=pltpar$mfrow)
for (i in seq_along(site_names)) {
  metrics = c('mete.50', 'Exp.50', 'Metric.50')
  true_fix = as.character(fixed$site) == site_names[i]
  true_log = as.character(logser$site) == site_names[i]
  index = match(site_names[i], shrtnames)
  true_fix = true_fix & fixed$area == area_of_interest[index]
  true_log = true_log & logser$area == area_of_interest[index]
  xlim = (range(c(0,fixed$Dist[true_fix]), na.rm=T))
  if (sum(true_log) > 0)
    ylim = (range(fixed[true_fix , metrics], logser[true_log, metrics], na.rm=T))
  else
    ylim = (range(fixed[true_fix , metrics], na.rm=T))
  x = ifelse(ylim[1] < 0, -1, 1)
  ylim = c(floor((ylim[1] %% x) * 10) / 10, 
           ceiling((ylim[2] %% 1) * 10) / 10)
  plot(Metric.50 ~ Dist, data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='o', lwd=5, cex=2, pch=19,
       xlab='', ylab='', frame.plot=F, axes=F, log='')
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0)
  hab_type = habitat[match(site_names[i], shrtnm)]
  mtext(side=3, paste(site_titles[i], '-', hab_type), cex=2)
  mtext(side=3, paste('(', LETTERS[i], ')', sep=''), adj=0, cex=2, font=2)
  metrics = metrics[-3]
  for (j in seq_along(metrics)) {
    lines(fixed[true_fix, 'Dist'], fixed[true_fix, metrics[j]],
          lwd=5, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 13)
    legend('bottomleft', 
           c('observed', 'recursive, observed SAD',
             'recursive, METE SAD','random, observed SAD'),
           col=c('black', rep(col, each=2)), cex=3, bty='n',
           lwd=8, lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()

png('./figs/ddr_arithY_logX_by_sites.png',
    width=480 * pltpar$wd_mult, height=480 * pltpar$ht_mult)
par(mfrow=pltpar$mfrow)
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
       xlim=xlim, ylim=ylim, type='o', lwd=5, cex=2, pch=19,
       xlab='', ylab='', frame.plot=F, axes=F, log='x')
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
          lwd=5, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 13)
    legend('bottomleft', 
           c('observed', 'recursive, observed SAD',
             'recursive, METE SAD','random, observed SAD'),
           col=c('black', rep(col, each=2)), cex=3, bty='n',
           lwd=rep(8, 5), lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()


png('./figs/ddr_loglog_by_sites.png',
    width=480 * pltpar$wd_mult, height=480 * pltpar$ht_mult)
par(mfrow=pltpar$mfrow)
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
       xlim=xlim, ylim=ylim, type='o', lwd=5, cex=2, pch=19,
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
          lwd=5, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 13)
    legend('bottomleft', 
            c('observed', 'recursive, observed SAD',
              'recursive, METE SAD','random, observed SAD'),
            col=c('black', rep(col, each=2)), cex=3, bty='n',
            lwd=rep(8, 5), lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()

png('./figs/ddr_loglog_by_sites_med.png',
    width=480 * pltpar$wd_mult, height=480 * pltpar$ht_mult)
par(mfrow=pltpar$mfrow)
for (i in seq_along(site_names)) {
  metrics = c('mete.50', 'Exp.50', 'Metric.50')
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
  plot(Metric.50 ~ Dist, data=fixed[true_fix, ],
       xlim=xlim, ylim=ylim, type='o', lwd=5, cex=2, pch=19,
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
          lwd=5, lty=lty[1], col=col[j], type='l')
    if (j == 1 & sum(true_log) > 0) { ## logser results only for METE model
      lines(logser[true_log, 'Dist'], logser[true_log, metrics[j]],
            lwd=5, lty=lty[2], col=col[j], type='l')
    }  
  }
  if(i == 13)
    legend('bottomleft', 
           c('observed', 'recursive, observed SAD',
             'recursive, METE SAD','random, observed SAD'),
           col=c('black', rep(col, each=2)), cex=3, bty='n',
           lwd=rep(8, 5), lty=c(1, lty), pch=c(19, rep(NA, 4)))
}
dev.off()

## Figure 3: one-to-one plots---------------------------------------------------
titles = c('recursive, METE SAD', 'recursive, observed SAD',
           'random, observed SAD')

cutoff = 10

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
    true = logser$indiv >= cutoff
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  else {
    true = fixed$indiv >= cutoff
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

png('./figs/ddr_one_to_one_median.png',
    width = 480 * 3, height= 480 * 1)
par(mfrow=c(1,3))
lims = c(0, .85)
for(i in 1:3) {
  metrics = c('mete.50', 'mete.50', 'Exp.50')
  if (i == 1) {
    x = logser[ , metrics[i]]
    y = logser[ , 'Metric.50']
  }
  else {
    x = fixed[ , metrics[i]]
    y = fixed[ , 'Metric.50']
  }
  plot((x), (y), type='n', axes=F, frame.plot=F, xlab='', ylab='',
       xlim=lims, ylim=lims)
  if (i == 1) {
    true = logser$indiv >= cutoff
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  else {
    true = fixed$indiv >= cutoff
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
    true = logser$indiv >= cutoff
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  else {
    true = fixed$indiv >= cutoff
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

png('./figs/ddr_loglog_one_to_one_median.png',
    width = 480 * 3, height= 480 * 1)
par(mfrow=c(1,3))
lims = c(2^-7, 1)
for(i in 1:3) {
  metrics = c('mete.50', 'mete.50', 'Exp.50')
  if (i == 1) {
    x = logser[ , metrics[i]]
    y = logser[ , 'Metric.50']
  }
  else {
    x = fixed[ , metrics[i]]
    y = fixed[ , 'Metric.50']
  }
  plot((x), (y), type='n', axes=F, frame.plot=F, xlab='', ylab='',
       xlim=lims, ylim=lims, log='xy')
  if (i == 1) {
    true = logser$indiv >= cutoff
    points((x[true]), (y[true]), pch=19)
    points((x[!true]), (y[!true]), pch=19, col='grey')
  } 
  else {
    true = fixed$indiv >= cutoff
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

## Supplemental site specific empirical DDR plots --------------------------
## generate pdfs
load('./sorensen/empirSorBin_bisect.Rdata')
load('./sorensen/empirSorAbu_bisect.Rdata')
load('simulated_empirical_results_bisect.Rdata')

for (metric in c('Sor')) {
  for (type in c('Bin', 'Abu')) {
    for (log_it in c(FALSE, TRUE)) {
      if (log_it)
        prefix = './figs/spat_empir_loglog_bisect_'
      else
        prefix = './figs/spat_empir_bisect_'
      pdf(paste(prefix, metric, '_', type, '_curves.pdf', sep=''),
          width = 7 * 3, height = 7 * 1)
      par(mfrow=c(1, 2))
      data = eval(parse(text=paste('empir', metric, type, sep='')))
      if (log_it) {
        plotEmpir(data, quants=FALSE, type='o', log='xy')
        plotEmpir(data, quants=TRUE, type='o', log='xy')
      }  
      else {
        plotEmpir(data, quants=FALSE, type='o')
        plotEmpir(data, quants=TRUE, type='o')
      }
      rm(data)
      dev.off()
    }  
  }  
}

## Supplemental site specific residual DDR plots -----------------------------------
load('./sorensen/empirSorAbu_bisect.Rdata')
load('simulated_empirical_results_bisect.Rdata')

resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)
resSorAbuLogSer = get_ddr_resid(empirSorAbu, simSorAbuLogSer)

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


## log-arith residual plot with lowess lines to summarize results across grains
png('./figs/ddr_arith_resid_by_sites_lowess_across_scales.png', width=480 * 4, height=480 * 4)
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
       xlab='', ylab='', frame.plot=F, axes=F)
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

## Sup Fig - SAR differences-----------------------------------------------
## Examine the difference in the SAR prediction between the recursive
## and semi-recursive formulations of METE

## NOTE: to run this code you must have first run './sar_run_all.R'
## load data

mk_sup_figs = FALSE
if (mk_sup_figs) {
  source('./scripts/spat_sar_load_and_avg_data.R')
  
  site_names = as.character(as.matrix(read.table('./data/shrtnames.txt'))) 
  site_names = "bci, sherman1, cocoli1, luquillo, bryan, bigoak, oosting, rocky, bormann, woodbridge, baldmnt, landsend, graveyard, ferp, serp, cross"
  site_names = unlist(strsplit(site_names, split=', '))
  site_titles = sub('1', '', site_names)
  site_titles = capwords(site_titles)
  site_titles[1] = "BCI"
  site_titles[11] = "Bald Mtn."
  site_titles[14] = "UCSC"
  site_titles[15] = "Serpentine"
  site_titles[16] = "Cross Timbers"
  
  
  png('./figs/sup_fig_mete_sim_analy_logser_sar_predictions.png',
      width=480 * 2, height=480 * 2)
  par(mfrow=c(4,4))
  for (i in seq_along(site_names)) {
    sar_tmp = subset(sar_data, subset=site==site_names[i])
    ylim = range(sar_tmp$logser_iter,
                 sar_tmp$logser_noniter,
                 sar_tmp$logser_avg,
                 na.rm=T)
    plot(logser_iter ~ area, data=sar_tmp, log='xy', 
         ylim=ylim,  type='n', main=site_titles[i],
         xlab='Area (m2)', ylab='Richness')
    ## Simulated log series
    if (site_names[i] != 'cross') {
      dat = meteAvgLogSer[[match(site_names[i], names(meteAvgLogSer))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(logser_avg ~ area, data=sar_tmp, lwd=3, col=1, lty=1)
    }
    ## Analytical Log series iterative
    lines(logser_iter ~ area, data=sar_tmp, lwd=3, col='red3')
    ## Analytical Log series noniterative
    lines(logser_noniter ~ area, data=sar_tmp, lwd=3, col='lightskyblue3')  
    ## Observed richness
    points(richness ~ area, data=sar_tmp, pch=1, cex=1.25)
    if (i == 1) {
      txt = c('Observed', 'Semi-recurisve CI', 'Semi-recursive Exp.',
              'Recursive Exp.', 'Non-recursive Exp.')
      legend('bottomright', txt, col=c('black', 'grey','black','red3', 'lightskyblue3'),
             lty=c(NA, rep(1,4)), lwd=c(NA, 6, rep(3,3)), pch=c(1, rep(NA,4)),
             bty='n', cex=1.25)
    }  
  }
  dev.off()
  
  ## arithmetic SARs
  png('./figs/sup_fig_mete_sim_analy_logser_sar_predictions_arith.png',
      width=480 * 2, height=480 * 2)
  par(mfrow=c(4,4))
  for (i in seq_along(site_names)) {
    sar_tmp = subset(sar_data, subset=site==site_names[i])
    ylim = range(sar_tmp$logser_iter,
                 sar_tmp$logser_noniter,
                 sar_tmp$logser_avg,
                 na.rm=T)
    plot(logser_iter ~ area, data=sar_tmp, log='', 
         ylim=ylim,  type='n', main=site_titles[i],
         xlab='Area (m2)', ylab='Richness')
    ## Simulated log series
    if (site_names[i] != 'cross') {
      dat = meteAvgLogSer[[match(site_names[i], names(meteAvgLogSer))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(logser_avg ~ area, data=sar_tmp, lwd=3, col=1, lty=1)
    }
    ## Analytical Log series iterative
    lines(logser_iter ~ area, data=sar_tmp, lwd=3, col='red3')
    ## Analytical Log series noniterative
    lines(logser_noniter ~ area, data=sar_tmp, lwd=3, col='lightskyblue3')  
    ## Observed richness
    points(richness ~ area, data=sar_tmp, pch=1, cex=1.25)
    if (i == 1) {
      txt = c('Observed', 'Semi-recurisve CI', 'Semi-recursive Exp.',
              'Recursive Exp.', 'Non-recursive Exp.')
      legend('bottomright', txt, col=c('black', 'grey','black','red3', 'lightskyblue3'),
             lty=c(NA, rep(1,4)), lwd=c(NA, 6, rep(3,3)), pch=c(1, rep(NA,4)),
             bty='n', cex=1.25)
    }  
  }
  dev.off()
  
}
## Sup Fig - Power vs Exponential models of DDR-------------------------------
## Compare the fit bettween the power and exponential models
## for both the METE predicted DDR and the observed DDR
source('./scripts/spat_functions.R')

load('./sorensen/empirSorAbu_bisect.Rdata')
load('simulated_empirical_results_bisect.Rdata')

## fix site names and ordering
index = sapply(names(empirSorAbu), 
               function(x) grep(x, names(simSorAbuFixed)))
simSorAbuFixed = simSorAbuFixed[index]
names(simSorAbuFixed) = names(empirSorAbu)

stats = list(mete = getStats(simSorAbuFixed, 'average'),
             empir = getStats(empirSorAbu, 'average'))

r2pwr = r2exp  = list(mete = NULL, empir = NULL)
for(i in seq_along(stats)) {
  r2pwr[[i]] = unlist(sapply(stats[[i]], function(x) x['pwr', 'r2', 'wtr',]))
  r2exp[[i]] = unlist(sapply(stats[[i]], function(x) x['exp', 'r2', 'wtr',]))
}

b1pwr = b0pwr = area = sites = list(mete = NULL, empir = NULL)
for(i in seq_along(stats)) {
  b1pwr[[i]] = unlist(sapply(stats[[i]], function(x) x['pwr', 'b1', 'wtr',]))
  b0pwr[[i]] = unlist(sapply(stats[[i]], function(x) x['exp', 'b0', 'wtr',]))
  area[[i]] = as.numeric(unlist(sapply(stats[[i]],
                                       function(x) names(x['pwr', 'b1', 'wtr',]))))
  n = unlist(lapply(stats[[i]], function(x)
                                 length(x['pwr', 'b1', 'wtr',])))
  sites[[i]] = unlist(mapply(rep, names(stats[[i]]), each=n))
}

## compute density kernals
dpwr = lapply(r2pwr, density)
dexp = lapply(r2exp, density)

png('./figs/sup_fig_r2_density_kernals.png', width=480*2, height=480)
par(mfrow=c(1,2))
linelwd = 3
xlims = list(c(.8, 1.01), c(0, 1.01))
ylims = list(c(0, 1250), c(0, 10))
title = c('METE DDR Functional Form', 'Empirical DDR Functional Form')
for(i in seq_along(dpwr)) {
  plot(dpwr[[i]], xlim=xlims[[i]], ylim=ylims[[i]], lwd=linelwd, col='black',
       main='',  xlab='', ylab='', frame.plot=F, axes=F)
  lines(dexp[[i]],  lwd=linelwd, col='grey')
  mtext(side=3, title[i], cex=2)
  addAxes()
  addylab('Density', padj=-2)
  addxlab(expression('Coefficient of Determination, ' * italic(R^2)),
          padj=2)
  if(i == 1)
    legend('topleft', c('Power', 'Exponential'), col=c('black','grey'),
           lwd = 5, lty=1, bty='n', cex=2)
} 
dev.off()

## set up graphic parameters 
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
shrtnm = shrtnm[-c(3, 6, 7)]
habitat = habitat[-c(3, 6, 7)]
hab = c('tropical', 'oak-hickory', 'pine', 'oak woodland', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
           "#FF8000")
col = habcol[match(habitat, hab)]

png('./figs/supl_fig_pwr_model_coef_vs_area.png', width=480*2, height=480*2)
par(mfrow=c(2,2))
title = c('METE DDR Functional Form', 'Empirical DDR Functional Form')
for(i in seq_along(area)) {
  plot(area[[i]], b0pwr[[i]], log='x', type='n', xlab='', ylab='',
       ylim=range(range(b0pwr), 0), axes=F)
  uni_sites = unique(sites[[i]])
  abline(h=0, col='grey', lwd=2, lty=2)
  for(j in seq_along(uni_sites)) {
    lines(area[[i]][sites[[i]] == uni_sites[j]],
          b0pwr[[i]][sites[[i]] == uni_sites[j]],
          col=col[j], lwd=3)
  }
  mtext(side=3, title[i], cex=2)
  addAxes()
  addylab(expression('Power model y-intercept, '* italic(beta[0])), padj=-.9)
  addxlab(expression('Area, ' * m^2), padj=1.75)
  if (i == 2)
    legend('bottomright', hab, lwd=6, col=habcol, lty=1, bty='n', cex=2)
}
for(i in seq_along(area)) {
  plot(area[[i]], b1pwr[[i]], log='x', type='n', xlab='', ylab='', 
       ylim=range(range(b1pwr), 0), axes=F)
  uni_sites = unique(sites[[i]])
  abline(h=0, col='grey', lwd=2, lty=2)
  for(j in seq_along(uni_sites)) {
    lines(area[[i]][sites[[i]] == uni_sites[j]],
          b1pwr[[i]][sites[[i]] == uni_sites[j]],
          col=col[j], lwd=3)
  }
  addAxes()
  addylab(expression('Power model slope, '* italic(beta[1])), padj=-.9)
  addxlab(expression('Area, ' * m^2), padj=1.75)
}
dev.off()


## Sup Fig - DDR Residuals Habitat Comparison----------------------------------

load('./sorensen/empirSorAbu_bisect.Rdata')
load('simulated_empirical_results_bisect.Rdata')

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

## bring in habitat type
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
dat$hab = habitat[match(dat$site, shrtnm)]
sites = unique(dat$site)

png('./figs/sup_fig_raw_ddr_residuals.png', width=480 * 2, height=480)
## raw emprical DDR residuals
par(mfrow=c(1,2))
for(j in 1:2){
  if(j == 1)
    main = 'METE'
  else
    main = 'RP'
  plot(avg.res ~ Dist, data=dat, log='x', type='n', ylim = c(-.5, .5), xlim=c(.5,512),
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

# Sup Fig - Analytical & Simulation DDR Comparison---------------------------
## this be take a few minutes to run
## requires compiled code that works on pc or ubuntu 

source('./scripts/spat_functions.R')
source('./scripts/spat_analytical_prob_funcs.R')
load_heap('./scripts')
dir.create('./tests')

## set parameters of METE community simulation
S = 100        ## number of species
n0 = 50      ## abundance of each species
N = n0 * S    ## total abundance
ncomm = 200   ## number of communites to generate
bisect = 8    ## number of bisections
shrtname = 'test'

sad = matrix(rep(n0, S), ncol=S)
write.table(sad, file='./tests/ddr_sad.csv', sep=',',
            row.names=FALSE, col.names=FALSE)

cmd = paste('python ./scripts/spat_heap_ddr.py', bisect, bisect, 'empirSAD',
            './tests/ddr_sad.csv', './tests/ddr_sor.csv')
system(cmd, wait=T)
analy_ddr = read.csv('./tests/ddr_sor.csv')

areas = 2^(0:bisect)
analy_sar = sapply(areas, function(A) S * (1 - heap_prob(0, A, n0, 2^8, use_c=T)))

cmd = paste('python ./scripts/spat_community_generation.py', S, N, ncomm, bisect,
            'False ./tests/ddr_sad.csv', shrtname)
system(cmd, wait=T)
comms = read.csv(paste('./comms/simulated_comms_', shrtname, '_empirSAD_C', ncomm,
                       '_B', bisect, '_grid.txt', sep=''))
comms = as.matrix(comms)

siml_sar = matrix(NA, nrow=ncomm, ncol=bisect + 1)
siml_shared = matrix(NA, nrow=ncomm, ncol=bisect)
siml_ddr = matrix(NA, nrow=ncomm, ncol=bisect)
for(i in 1:ncomm) {
  tmp_comm = comms[comms[,1] == i, ]
  sp_mat = tmp_comm[ , -(1:3)]
  coords = tmp_comm[, 2:3]
  psp = mat2psp(sp_mat, coords)
  sar = getSAR(psp, grains=areas)
  shared = vario_bisect(sp_mat > 0,  coords, sep_orders=1:bisect,
                        distance.metric='shared')
  sor = vario_bisect(sp_mat > 0,  coords, sep_orders=1:bisect,
                     distance.metric='bray', NA_replace=0)
  ## export results
  siml_sar[i, ] = data.frame(sar)$richness
  siml_shared[i, ] = shared$vario$exp.var
  siml_ddr[i, ] = 1 - sor$vario$exp.var
}

avg_sar = apply(siml_sar, 2, mean)
avg_shared = apply(siml_shared, 2, mean)
avg_ddr = apply(siml_ddr, 2, mean)

qts_sar = apply(siml_sar, 2, quantile, c(.025, .975))
qts_shared = apply(siml_shared, 2, quantile, c(.025, .975))
qts_ddr = apply(siml_ddr, 2, quantile, c(.025, .975))

png('./figs/sup_fig_analy_simul_sar_ddr_check.png', width=480*3 ,
    height=480)
par(mfrow=c(1, 3))
## compare SARs
plot(areas, avg_sar, type='n', xlab='', ylab='', log='xy', axes=F,
     ylim=range(qts_sar))
addxlab('Area', padj=2)
addylab('Richness', padj=-1)
addCI(areas, qts_sar[1, ], qts_sar[2, ], col='pink')
lines(areas, avg_sar, lwd=3, col='red')
lines(areas, analy_sar, lwd=3, col='blue')
addAxes()
## compare # shared species versus distance
plot(analy_ddr$dist, avg_shared, type='n', xlab='', ylab='',
     ylim=range(qts_shared), axes=F)
addxlab('Distance', padj=2)
addylab('Number of shared species', padj=-1)
addCI(analy_ddr$dist, qts_shared[1, ], qts_shared[2, ], col='pink')
lines(analy_ddr$dist, avg_shared, col='red', lwd=3)
lines(analy_ddr$dist, analy_ddr$sor * analy_sar[1], col='blue', lwd=3)
addAxes()
## compare DDR curves
plot(analy_ddr$dist, avg_ddr, type='n', ylim=range(qts_ddr),
     xlab='',ylab='', axes=F)
addxlab('Distance', padj=2)
addylab('Sorensen Similarity', padj=-1)
addCI(analy_ddr$dist, qts_ddr[1, ], qts_ddr[2, ], col='pink')
lines(analy_ddr$dist, avg_ddr, col='red', lwd=3)
lines(analy_ddr$dist, analy_ddr$sor, col='blue', lwd=3)
legend('topright',c('Simul. CI', 'Simul. Avg.','Analy. Exp.'),
       col=c('pink', 'red','blue'), lty=1, lwd=c(10, 5, 5),
       bty='n', cex=3)
addAxes()

dev.off()

print('Generating figures, complete!')

