## METE analytical SAR predictions vs Empirical SARs

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

## load data
source('./scripts/spat_sar_load_and_avg_data.R')

## plot an example SAR figure for presentation ---------------------------------
par(mfrow=c(1,1))
site = 'bigoak'
purple = rgb(112, 48, 160, maxColorValue=160)
lightblue = "#1AB2FF"

i = match(site, names(meteEmpirSAD))
  ## log-log
    plot(sr_iter ~ area, data=meteEmpirSAD[[i]], 
         ylim=range(c(meteEmpirSAD[[i]]$sr_iter, empir[[i]]$richness)),
         xlim=range(c(meteEmpirSAD[[i]]$area, empir[[i]]$area)), log='xy',
         type='n', main=names(meteEmpirSAD)[i], ylab='SR', xlab='Area (m2)')
    ## meteEmpirSAD CI
    dat = meteAvgEmpirSAD[[match(names(meteEmpirSAD)[i], names(meteAvgEmpirSAD))]]
#    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
    lines(sr.avg ~ grains, data=dat, lwd=3, col='grey', lty=2)
    ## RP CI
    dat = srExp[[match(names(meteEmpirSAD)[i], names(srExp))]]
#    addCI('grains', 'S_lo', 'S_hi', data='dat', col='pink')
    lines(S_binom ~ grains, data=dat, col='red', lwd=3)
    ## analytical meteEmpirSAD    
    lines(sr_noniter ~ area, data=meteEmpirSAD[[i]], col='dodgerblue', lwd=3)
    lines(sr_iter ~ area, data=meteEmpirSAD[[i]], col='grey', lwd=3, lty=1)
    ## data
    lines(richness ~ area, data = empir[[i]], pch=19, type='o', lwd=3)
    legend('bottomright', c('Empirical','RP','meteEmpirSAD sim', 'meteEmpirSAD noniter'),
           pch=c(19, rep(NA, 3)), col=c(1, 'red', 'grey', 'dodgerblue'),
           bty='n', lwd=3)



## Compare METE SAR predictions--------------------------------------------
pdf('./figs/mete_sar_predictions.pdf', width=7 * 2, height=7 * 2)
  par(mfrow=c(4,4))
  for (i in seq_along(meteLogSer)) {
    plot(sr_iter ~ area, data=meteLogSer[[i]], log='xy',
         ylim=range(c(meteLogSer[[i]]$sr_iter, empir[[i]]$richness)),
         xlim=range(c(meteLogSer[[i]]$area, empir[[i]]$area)),
         type='n', main=names(meteLogSer)[i], xlab='Area (m2)')
    ## Simulated log series
    if (names(meteLogSer)[i] != 'cross') {
      dat = meteAvgLogSer[[match(names(meteLogSer)[i], names(meteAvgLogSer))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
    }
    ## Simulated fixed SAD
    dat = meteAvgEmpirSAD[[match(names(meteLogSer)[i], names(meteAvgEmpirSAD))]]
    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='green3')
    ## Analytical Log series iterative
    lines(sr_iter ~ area, data=meteLogSer[[i]], type='o', col='red')
    ## Analytical log series noniterative
    if (names(meteLogSer)[i] != 'cross')
      lines(sr_noniter ~ area, data=meteLogSer[[i]], type='o', pch=19, col='blue')
    ## Analytical fixed SAD iterative
    lines(sr_iter ~ area, data=meteEmpirSAD[[i]], type='o', col='red', lty=2)
    ## Analytical fixed SAD noniterative
    lines(sr_noniter ~ area, data=meteEmpirSAD[[i]], type='o', pch=19, col='blue', lty=2)
    if (i == 1) {
      txt = c('Sim CI LogSer', 'Sim CI EmpirSAD', 'Analy iter LogSer', 
              'Analy iter Empir SAD', 'Analy noniter LogSer', 'Analy noniter EmpirSAD')
      legend('bottomright', txt, col=c('grey','green3','red','red', 'blue', 'blue'),
             pch=c(NA, NA,1,1,19,19), lty=c(1,1,1,2,1,2), lwd=c(2,2, rep(1, 4)),
             bty='n', cex=1)    
    }  
  }
dev.off()

## Compare SAR models and data--------------------------------------------------

## for LogSer SAD
pdf('./figs/meteLogSer_&_empir_sar.pdf', width=7 * 2, height=7 * 2)
  par(mfrow=c(4,4))
  ## log-log
  for (i in seq_along(meteLogSer)) {
    plot(sr_iter ~ area, data=meteLogSer[[i]], 
         ylim=range(c(meteLogSer[[i]]$sr_iter, empir[[i]]$richness)),
         xlim=range(c(meteLogSer[[i]]$area, empir[[i]]$area)), log='xy',
         type='n', main=names(meteLogSer)[i], xlab='Area (m2)')
    ## meteLogSer CI
    if (names(meteLogSer)[i] != 'cross') {
      dat = meteAvgLogSer[[match(names(meteLogSer)[i], names(meteAvgLogSer))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(sr.avg ~ grains, data=dat, lwd=2, col='darkgrey')
    }  
    ## RP CI
    dat = srExp[[match(names(meteLogSer)[i], names(srExp))]]
    lines(S_logser_binom ~ grains, data=dat, col='red', lwd=2)
    ## analytical meteLogSer    
    lines(sr_iter ~ area, data=meteLogSer[[i]], type='o', col='green3')
    if (names(meteLogSer)[i] != 'cross') 
      lines(sr_noniter ~ area, data=meteLogSer[[i]], type='o', col='dodgerblue')
    ## data
    lines(richness ~ area, data = empir[[i]], pch=19, type='o')
    if(i == 1)
      legend('bottomright', c('Empirical','RP','meteLogSer sim', 'meteLogSer iter', 'meteLogSer noniter'),
             pch=c(19, NA, NA, 1, 1), col=c(1, 'red', 'grey', 'green3', 'dodgerblue'),
             bty='n', lwd=c(2, 3, 3, 2, 2))
  }
dev.off()    

## for EmpirSAD
pdf('./figs/meteEmpirSAD_&_empir_sar.pdf', width=7 * 2, height=7 * 2)
  par(mfrow=c(4,4))
  ## log-log
  for (i in seq_along(meteEmpirSAD)) {
    plot(sr_iter ~ area, data=meteEmpirSAD[[i]], 
         ylim=range(c(meteEmpirSAD[[i]]$sr_iter, empir[[i]]$richness)),
         xlim=range(c(meteEmpirSAD[[i]]$area, empir[[i]]$area)), log='xy',
         type='n', main=names(meteEmpirSAD)[i], xlab='Area (m2)')
    ## meteEmpirSAD CI
    dat = meteAvgEmpirSAD[[match(names(meteEmpirSAD)[i], names(meteAvgEmpirSAD))]]
    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
    lines(sr.avg ~ grains, data=dat, lwd=2, col='darkgrey')
    ## RP CI
    dat = srExp[[match(names(meteEmpirSAD)[i], names(srExp))]]
#    addCI('grains', 'S_lo', 'S_hi', data='dat', col='pink')
    lines(S_binom ~ grains, data=dat, col='red', lwd=2)
    ## analytical meteEmpirSAD    
    lines(sr_iter ~ area, data=meteEmpirSAD[[i]], type='o', col='green3')
    lines(sr_noniter ~ area, data=meteEmpirSAD[[i]], type='o', col='dodgerblue')
    ## data
    lines(richness ~ area, data = empir[[i]], pch=19, type='o')
    if(i == 1)
      legend('bottomright', c('Empirical','RP','meteEmpirSAD sim', 'meteEmpirSAD iter', 'meteEmpirSAD noniter'),
             pch=c(19, NA, NA, 1, 1), col=c(1, 'red', 'grey', 'green3', 'dodgerblue'),
             bty='n', lwd=c(2, 3, 3, 2, 2))
  }
dev.off()   

## compare residuals between METE and RP----------------------------------------

shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
col = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
        "#FF8000")

as.matrix(sort(apply(abs(sar_res[ , 3:10]), 2, mean, na.rm=T)))
as.matrix(sort(apply(abs(sar_res[!is.na(sar_res[,4]), 3:10]), 2, mean, na.rm=T)))

## average residuals on a per site basis
## then examine which model performs the best for each site
## this keeps sites with lots of spatial scales from dominating
## the residuals

avg_sar_res = aggregate(sar_res[ , -c(1:2, 11:12)] / sar_res$richness, by=list(sar_res$site), 
                        FUN = function(x) mean(x^2, na.rm=T))
indices = apply(avg_sar_res[ , -1], 1, function(x) which(min(x, na.rm=T) == x))
wins = as.matrix(table(names(avg_sar_res[ , -1])[c(indices,1:8)]) - 1)
res_avg = apply(avg_sar_res[ , -1], 2, mean, na.rm=T)[order(names(avg_sar_res[,-1]))]
cbind(wins, res_avg)
                    res_avg
empirsad_avg     0 15.75582
empirsad_iter    0 54.45989
empirsad_noniter 6  3.40553
empirsad_rp      1 28.08008
logser_avg       3 44.68045
logser_iter      2 59.86338
logser_noniter   2 27.78716
logser_rp        2 13.26835

## compare only at the empirical SAD models
avg_sar_res = aggregate(sar_res[ , 7:10] / sar_res$richness, by=list(sar_res$site), 
                        FUN = function(x) mean(x^2, na.rm=T))
indices = apply(avg_sar_res[ , -1], 1, function(x) which(min(x) == x))
wins = as.matrix(table(names(avg_sar_res[ , -1])[c(indices,1:4)]) - 1)
res_avg = apply(avg_sar_res[ , -1], 2, mean)[order(names(avg_sar_res[,-1]))]
cbind(wins, res_avg)
                     res_avg
empirsad_avg      3 15.75582
empirsad_iter     1 54.45989
empirsad_noniter 10  3.40553
empirsad_rp       2 28.08008

## normalized by Savg
                        res_avg
empirsad_avg      0 0.041606515
empirsad_iter     1 0.066459792
empirsad_noniter 14 0.007523939
empirsad_rp       1 0.032704148


## compare only at the empirical SAD models simulated iterative and rp
avg_sar_res = aggregate(sar_res[ , c('empirsad_avg','empirsad_rp')] ,
                        by=list(sar_res$site), 
                        FUN = function(x) mean(x^2, na.rm=T))
indices = apply(avg_sar_res[ , -1], 1, function(x) which(min(x) == x))
wins = as.matrix(table(names(avg_sar_res[ , -1])[c(indices,1:4)]) - 1)
res_avg = apply(avg_sar_res[ , -1], 2, mean)[order(names(avg_sar_res[,-1]))]
cbind(wins, res_avg)
                 res_avg
empirsad_avg 12 15.75582
empirsad_rp   4 28.08008

## normalized by Savg
                   res_avg
empirsad_avg  4 0.04160652
empirsad_rp  12 0.03270415


## compare empirical analytical iterative and rp
avg_sar_res = aggregate(sar_res[ , c('empirsad_iter','empirsad_rp')],
                        by=list(sar_res$site), 
                        FUN = function(x) mean(x^2, na.rm=T))
indices = apply(avg_sar_res[ , -1], 1, function(x) which(min(x) == x))
wins = as.matrix(table(names(avg_sar_res[ , -1])[c(indices,1:4)]) - 1)
res_avg = apply(avg_sar_res[ , -1], 2, mean)[order(names(avg_sar_res[,-1]))]
cbind(wins, res_avg)
                 res_avg
empirsad_iter  6 54.45989
empirsad_rp   10 28.08008


## plot residuals---------------------------------------------------------------
pdf('./figs/sar_residuals.pdf', width=7 * 2, height=7 * 2)
  sites = unique(sar_res$site)
  par(mfrow=c(2,2))
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE noniter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_noniter ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  legend('bottomright', hab, col=col, lwd=2, bty='n')
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_iter ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter sim')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_avg ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='RP')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_rp ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex], lty=1)
  }
dev.off()

## plot normalized residuals----------------------------------------------------
pdf('./figs/sar_norm_residuals.pdf', width=7 * 2, height=7 * 2)
  sites = unique(sar_res$site)
  par(mfrow=c(2,2))
  ylims = c(-.6, .6)  
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims,
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE noniter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_noniter / richness ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  legend('bottomright', hab, col=col, lwd=2, bty='n')
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims,
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_iter / richness ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims, 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter sim')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_avg / richness ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims, 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='RP')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_rp / richness ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex], lty=1)
  }
dev.off()

## plot log residuals----------------------------------------------------
pdf('./figs/sar_log_diff_residuals.pdf', width=7 * 2, height=7 * 2)
  sites = unique(sar_res$site)
  par(mfrow=c(2,2))
  ylims = c(-.7, .7)  
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims,
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE noniter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(log(richness / (-empirsad_noniter + richness)) ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  legend('bottomright', hab, col=col, lwd=2, bty='n')
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims,
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(log(richness / (-empirsad_iter + richness)) ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims, 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter sim')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(log(richness / (-empirsad_avg + richness)) ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims, 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='RP')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(log(richness / (-empirsad_rp + richness)) ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex], lty=1)
  }
dev.off()

## plot studentized residuals----------------------------------------------------
pdf('./figs/sar_stud_residuals.pdf', width=7 * 2, height=7 * 2)
  sites = unique(sar_res$site)
  par(mfrow=c(2,2))
  ylims = c(-6, 6)  
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims,
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE noniter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_noniter / std ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  legend('bottomright', hab, col=col, lwd=2, bty='n')
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims,
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_iter / std ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims, 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='METE iter sim')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_avg / std ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex])
  }
  plot(empirsad_noniter ~ area, data=sar_res, log='x', ylim=ylims, 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), main='RP')
  axis(side=1, cex.axis=1.75, padj=.5, lwd=8, at=10 ^ (-1:6))
  axis(side=2, cex.axis=1.75, lwd=8)
  abline(h=0, lwd=5)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_rp / std ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=col[habindex], lty=1)
  }
dev.off()


## plot universal curve --------------------------------------------------------

slopes = read.csv('./sar/sar_slopes.csv')
slopes = data.frame(slopes,logNS = log(slopes$NS))
## drop crosstimbers mete prediction
slopes$predZ[slopes$comm == 'cross'] = NA
## average cocoli and sherman
cocoliAvg = (slopes[slopes$comm == 'cocoli1', -1] + 
             slopes[slopes$comm == 'cocoli2', -1] ) / 2
shermanAvg = (slopes[slopes$comm == 'sherman1', -1] + 
              slopes[slopes$comm == 'sherman2', -1] ) / 2
slopes[slopes$comm == 'cocoli1', -1] = cocoliAvg
slopes[slopes$comm == 'sherman1', -1] = shermanAvg
slopes = slopes[!(slopes$comm %in% c('cocoli2','sherman2','sherman3')), ]

## attempt to censor the predicted zvals
## first sort with respect to log(N/S)
ratio = slopes$logNS
z = slopes$predZ[order(ratio)]
ratio = sort(ratio)
## drop NAvals
ratio = ratio[!is.na(z)]
z = z[!is.na(z)]
## drop z greater than one
ratio = ratio[!(z > 1)]
z = z[!(z > 1)]
## compute a smoother
plot(ratio, z)
Zsmooth = lowess(ratio,z, f=.2)
lines(Zsmooth, col='red', lwd=2)
res = abs(Zsmooth$y - z)

pdf('./figs/universal_sar.pdf', width=7 * 1.25, height=7)
  par(mfrow=c(1,1))
  plot(predZ ~ logNS, data=slopes, type='n', ylim=c(0,1), frame.plot=F,
       xlab='ln(N/S)', ylab='Z-value', xlim = c(0, 8))
  lines(Zsmooth, lwd=4, lty=2)
  comms = unique(slopes$comm)
  col = colorRampPalette(c("blue", "red", "green", "orange"))(length(comms))
  for (i in seq_along(comms))
    points(obsZ ~ logNS, data=slopes, subset=comm == comms[i], col=col[i], 
           pch=19, cex=1.25)
  legend('topright', as.character(comms), col=col, pch=19, bty='n', cex=1)
dev.off()

## for slide
pch = c(17, 0, 16, 1, 15, 2)
pch = rep(19, 6)

par(mfrow=c(1,1))
plot(predZ ~ logNS, data=slopes, type='n', ylim=c(0,1), axes=F, frame.plot=F,
     xlab='', ylab='')
axis(side=1,lwd=4,cex.axis=2)
axis(side=2,lwd=4,cex.axis=2)
lines(Zsmooth, lwd=4, lty=2)
comms = unique(slopes$comm)
for (i in seq_along(comms)) {
  habindex = match(habitat[match(comms[i], shrtnm)], hab)
  points(obsZ ~ logNS, data=slopes, subset=comm == comms[i], col=col[habindex], 
         pch=pch[habindex], cex=1.5, lwd=2)
}  

plot(1:10, 1:10, type='n', xlab='', ylab='', axes=F, frame.plot=F)
legend('center', hab, col=col, pch=pch, lwd=6, lty=NA, cex=3, bty='n')

## compute z-values for analytical iterative curve
sar = read.csv('./sar/param_sar_avgs.csv')
sar$z = NA
for (i in 1:nrow(sar)){
  if (sar$grains[i] < 256)
    sar$z[i] = log(sar$sr.avg[i+1] / sar$sr.avg[i]) / log(4)
}

pdf('./figs/param_sar_scale_collapse.pdf', width=7 * 2, height=7)
 par(mfrow=c(1,2))
 plot(log(sar$N/sar$S) , sar$z)
 lines(Zsmooth)
 plot(log(sar$N/sar$S) / log(4096 / sar$grains), sar$z)
 ## x-axis divided by two here b/c 2 bisections rel to anchor scale
 lines(Zsmooth$x / 2, Zsmooth$y) 
dev.off()
