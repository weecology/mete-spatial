setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

## Figure 2 - example Simulated DDR with simulation results---------------------
#pdf('./figs/fig2_example_ddr_w_sim_results.pdf', width=7 * 3, height=7)
#windows(width=7 * 3, height=7)
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


dev.off()
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

#pdf('./figs/fig3_sar_ddr_empirical.pdf', width= 7 * 3, height= 7 * 2)
#windows(width= 7 * 3, height= 7 * 2)

par(mfrow=c(2, 3))

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

## panels B & C: empirical SAR residuals
  sites = unique(sar_res$site)
  plot(empirsad_avg ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6))
  addAxis1(at=10 ^ seq(-1, 5, 2))
  addAxis2()
  abline(h=0, lwd=4, lty=3)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_avg ~ area, data=sar_res, subset= site == sites[i],
          lwd=3, col=habcol[habindex])
  }
#  legend('bottomright', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')
  ##
  plot(empirsad_avg ~ area, data=sar_res, log='x', ylim=c(-25, 25), 
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6))
  addAxis1(at=10 ^ seq(-1, 5, 2))
  addAxis2()
  abline(h=0, lwd=3, lty=3)
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
sar_data = read.csv('./sar/empir_sars.csv')
sar_data$area = round(sar_data$area, 2)
dat = merge(dat, sar_data[ , c('site', 'area', 'richness', 'indiv')], all.x=TRUE)
## subset so that has at least 20 individuals
dat = subset(dat, indiv >= 20)

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
    abline(h=0, lwd=4, lty=3)
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


dev.off()

## Supplemental Figure 1--------------------------------------------------------
## r2 of model fits to simulated results
load('./sorensen/simSorAbuAvg.Rdata')
S = round(10^seq(log10(10), log10(100), length.out=20))
N = round(10^seq(log10(120), log10(5e5), length.out=20))
#stats = getSimStats(simSorAbuAvg, S, N)

#pdf('./figs/sup_fig1_r2_simulated_ddr.pdf', width = 7, height= 7)
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
addAxis1()
addAxis2()
addxlab(expression('Coefficient of Determination, ' * italic(R^2)),
        padj=2)
addylab('Freqency')
lines(xpwr, ypwr,  lwd=linelwd, col='black')
lines(xexp, yexp, lwd=linelwd, col='grey')

legend('topleft', c('Exponential Model', 'Power  Model'),
       lty=1, lwd=6, bty='n', cex=1.5, col=c( 'grey', 'black'))


dev.off()

## Supplemental Fig 2. METE DDR scale collapse for all grains-----------------------------------------
ddr = read.csv('./sorensen/param_ddr_wtr_pwr_stats.csv')
grains = unique(ddr$grains)

#pdf('./figs/sup_fig2_mete_ddr_scale_collapse_all_grains.pdf', width= 7 * 2, height= 7)
#windows(width= 7 * 2, height= 7)

  par(mfrow=c(1,2))
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  ddr$ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)
##
  plot(10^b0 ~ ratio , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(0, 1))
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(italic(N(A))/italic(S(A)))),
          padj=2.25)
  mtext(side=2, expression('Initial Similarity, ' * hat(beta)[0]), 
        cex=1.75, padj=-1)
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio, 10^tmp$b0), col=col[g], lwd=4)
  }
##
  plot(b1 ~ ratio , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F)
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(italic(N(A))/italic(S(A)))),
          padj=2.25)
  mtext(side=2, expression('Decay Rate, ' * hat(beta)[1]), 
        cex=1.75, padj=-1)
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    lines(lowess(tmp$ratio, tmp$b1), col=col[g], lwd=4)
  }  
  legend('topright', legend=c('Grains',  unique(ddr$grains)), cex=2,
         col=c(NA, col), lwd=c(NA, rep(4, 5)), bty='n')

dev.off()

## Supplemental Fig 3. METE DDR scale collapse with data---------------------------------------------
## read in stats for the simulated parameter space
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
    ratio = log(N0 / S0) / log(navg / savg)
    if (exists('dd_stats'))
      dd_stats = rbind(dd_stats, 
                       data.frame(site, area=grains[g], A0, S0, N0, savg, navg, b0, b1, ratio))
    else
      dd_stats = data.frame(site, area=grains[g], A0, S0, N0, savg, navg, b0, b1, ratio)
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
#pdf('./figs/sup_fig3_mete_ddr_scale_collapse_with_empir_data.pdf', width= 7 * 2, height=7 * 1)
#windows(width = 7 * 2, height=7 * 1)

ddr$ratio = log(ddr$N / ddr$S) / log(ddr$ind.avg / ddr$sr.avg)

par(mfrow=c(1,2))
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  plot(10^b0 ~ ratio, data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, log='', ylim=c(0, 2))
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(italic(N(A))/italic(S(A)))),
          padj=2.25)
  mtext(side=2, expression('Initial Similarity, ' * hat(beta)[0]), 
        cex=1.75, padj=-1)
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    points(tmp$ratio, 10^tmp$b0, col='grey', lwd=1, pch=19)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    points(10^b0 ~ ratio, data=dd_stats, subset= site == sites[i] & navg >= 20, 
           col=habcol[habindex], lwd=2, lty=2)
  }
##
  plot(b1 ~ ratio, data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(-.7, .1), log='')
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(italic(N(A))/italic(S(A)))),
          padj=2.25)
  mtext(side=2, expression('Decay Rate, ' * hat(beta)[1]), 
        cex=1.75, padj=-1) 
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    points(tmp$ratio, tmp$b1, col='grey', lwd=1, pch=19)
#    lines(lowess(tmp$ratio, tmp$b1), col=col[g], lwd=4)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    points(b1 ~ ratio, data=dd_stats, subset= site == sites[i] & navg >= 20, 
           col=habcol[habindex], lwd=2, lty=2)
  }
  legend('topright', c('METE Simulated', hab), cex=1.5,
         col=c('grey', habcol), lty=NA, lwd=3, pch=c(19, rep(1,5)), bty='n')


dev.off()

##Alternative x-axis for slope collapse--------------

#pdf('./figs/sup_fig3_mete_ddr_scale_collapse_with_empir_data_alt_scaling.pdf', width= 7 * 2, height=7 * 1)
#windows(width = 7 * 2, height=7 * 1)

par(mfrow=c(1,2))
  col = colorRampPalette(c('dodgerblue', 'red'))(5)
  plot(10^b0 ~ ratio , data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, log='', ylim=c(0, 2))
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(italic(N(A))/italic(S(A)))),
          padj=2.25)
  mtext(side=2, expression('Initial Similarity, ' * hat(beta)[0]), 
        cex=1.75, padj=-1)
  for (g in seq_along(grains)) { 
    tmp = subset(ddr, grains == grains[g])
    points(tmp$ratio, 10^tmp$b0, col='grey', lwd=1, pch=19)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    points(10^b0 ~ ratio, data=dd_stats, subset= site == sites[i] & navg >= 20, 
           col=habcol[habindex], lwd=2, lty=2)
  }
##
  plot(b1 ~ I(ratio / log(4096 / grains)), data=ddr, xlab='', ylab='', type='n',
       frame.plot=F, axes=F, ylim=c(-.7, .1), log='')
  addAxis1()
  addAxis2()
  addxlab(expression(log(italic(N)[0]/italic(S)[0]) / log(italic(N(A))/italic(S(A))) / log(italic(A[0])/italic(A))),
          padj=2.25)
  mtext(side=2, expression('Decay Rate, ' * hat(beta)[1]), 
        cex=1.75, padj=-1) 
  for (g in seq_along(grains)) {
    tmp = subset(ddr, grains == grains[g])
    points(tmp$ratio / log(4096 / tmp$grains), tmp$b1, col='grey', lwd=1, pch=19)
#    lines(lowess(tmp$ratio, tmp$b1), col=col[g], lwd=4)
  }  
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    points(b1 ~ I(ratio / log(A0 / A)), data=dd_stats, subset= site == sites[i] & navg >= 20, 
           col=habcol[habindex], lwd=2, lty=2)
  }
  legend('topright', c('METE Simulated', hab), cex=1.5,
         col=c('grey', habcol), lty=NA, lwd=3, pch=c(19, rep(1,5)), bty='n')


dev.off()

## Supplemental Fig 4. Normalized SAR Residuals-------------------------------------------

source('./scripts/spat_sar_load_and_avg_data.R')
## set up graphic parameters 
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
habcol = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
        "#FF8000")

#pdf('./figs/sup_fig4_normalized_sar_residuals.pdf', width=7 * 2, height=7)
#windows(width=7 * 2, height=7)

par(mfrow=c(1,2))
  sites = unique(sar_res$site)
  plot(empirsad_avg / richness ~ area, data=sar_res, log='x', 
       type='n', frame.plot=F, axes=F, xlab='', ylab='', 
       xlim = c(0.1, 1e6), ylim=c(-.6, .6))
  mtext(side=3, 'METE', cex=2)
  addAxis1(at=10 ^ seq(-1, 5, 2))
  addAxis2()
  addxlab(expression('Area ' * (m^2)), padj=2)
  addylab('Normalized Residuals')
  abline(h=0, lwd=4, lty=3)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_avg / richness~ area, data=sar_res, subset= site == sites[i],
          lwd=3, col=habcol[habindex])
  }
#  legend('bottomright', hab, col=habcol, lty=1, lwd=7, cex=2, bty='n')
  ##
  plot(empirsad_avg / richness ~ area, data=sar_res, log='x',
       type='n', frame.plot=F, axes=F, xlab='', ylab='',
       xlim = c(0.1, 1e6), ylim=c(-.6, .6))
  mtext(side=3, 'Random Placement', cex=2)
  addAxis1(at=10 ^ seq(-1, 5, 2))
  addAxis2()
  addxlab(expression('Area ' * (m^2)), padj=2)
  addylab('Normalized Residuals')
  abline(h=0, lwd=4, lty=3)
  for (i in seq_along(sites)) {
    habindex = match(habitat[match(sites[i], shrtnm)], hab)
    lines(empirsad_rp / richness ~ area, data=sar_res, subset= site == sites[i],
          lwd=4, col=habcol[habindex], lty=1)
  }
  legend('topright', hab, col=habcol, lty=1, lwd=5, bty='n', cex=1.25)

dev.off()
