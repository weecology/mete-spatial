setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

## load data
source('./scripts/spat_sar_load_and_avg_data.R')

indices = c(4:5, 8:9)

pred_sar = sar_res$richness - sar_res[ , indices]
pred_sar = pred_sar[ , c(1, 3, 2, 4)]
pred_sar = data.frame(pred_sar, richness = sar_res$richness, 
                      area = sar_res$area, site = sar_res$site,
                      hab = sar_res$hab)


## plot individual site relationships-------------------------------------------

site_names = "bci, sherman1, cocoli1, luquillo, bryan, bigoak, oosting, rocky, bormann, woodbridge, baldmnt, landsend, graveyard, ferp, serp, cross"
site_names = unlist(strsplit(site_names, split=', '))

col = c('red', 'red', 'dodgerblue', 'dodgerblue')
lty = c(1, 3, 1, 3)

## mete 
jpeg('./figs/mete_&_empir_sar.jpeg', width=480 * 4, height=480 * 4, quality=100)
#pdf('./figs/mete_&_empir_sar.pdf', width=7 * 4, height=7 * 4)
  par(mfrow=c(4,4))
  ## log-log
  for (i in seq_along(site_names)) {
    true = as.character(pred_sar$site) == site_names[i]
    xlim = range(pred_sar$area[true], na.rm=T)
    ylim = range(pred_sar[true , 1:5], na.rm=T)
    xlim = 2^c(floor(log2(xlim[1])), ceiling(log2(xlim[2])))
    ylim = 2^c(floor(log2(ylim[1])), ceiling(log2(ylim[2])))
    plot(richness ~ area, data=pred_sar[true, ],
         xlim=xlim, ylim=ylim, log='xy', type='n',
         xlab='', ylab='', frame.plot=F, axes=F)
    xticks = 2^seq(log2(xlim[1]), log2(xlim[2]), 2)
    yticks = 2^seq(log2(ylim[1]), log2(ylim[2]), 2)
    addAxis(side=1, cex.axis=3, padj=.75, at=xticks, lab=xticks)
    addAxis(side=2, cex.axis=3, padj=0, at=xticks, lab=xticks)
    mtext(side=3, unique(pred_sar$hab[true]), cex=2)
    mtext(side=3, paste(letters[i], ')', sep=''), adj=0, cex=3)
    for (j in 1:4) 
      lines(pred_sar[true, j] ~ area, data=pred_sar[true, ], col=col[j], lwd=4,
            lty=lty[j])
    ## data
    points(richness ~ area, data = pred_sar[true, ], pch=1, cex=3, lwd=2)
    if(i == 1)
      legend('bottomright', 
             c('observed', 'iterative', 'iterative observed SAD', 'noniterative',
               'noniterative observed SAD'), pch=c(1, rep(NA, 4)), col=c(1, col),
             cex=2.5, bty='n', lwd=c(2, rep(4, 4)), lty=c(NA, lty))
  }
dev.off()

jpeg('./figs/mete_&_empir_sar_log2.jpeg', width=480 * 4, height=480 * 4, quality=100)
#pdf('./figs/mete_&_empir_sar_log2.pdf', width=7 * 4, height=7 * 4)
par(mfrow=c(4,4))
## log2-log2
for (i in seq_along(site_names)) {
  true = as.character(pred_sar$site) == site_names[i]
  xlim = log2(range(pred_sar$area[true], na.rm=T))
  ylim = log2(range(pred_sar[true , 1:5], na.rm=T))
  xlim = c(floor(xlim[1]), ceiling(xlim[2]))
  ylim = c(floor(ylim[1]), ceiling(ylim[2]))
  plot(log2(richness) ~ log2(area), data=pred_sar[true, ],
       xlim=xlim, ylim=ylim, type='n',
       xlab='', ylab='', frame.plot=F, axes=F)
  addAxis(side=1, cex.axis=3, padj=.75)
  addAxis(side=2, cex.axis=3, padj=0)
  mtext(side=3, unique(pred_sar$hab[true]), cex=2)
  mtext(side=3, paste(letters[i], ')', sep=''), adj=0, cex=3)
  for (j in 1:4) 
    lines(log2(pred_sar[true, j]) ~ log2(area), data=pred_sar[true, ], col=col[j],
          lwd=4, lty=lty[j])
  ## data
  points(log2(richness) ~ log2(area), data = pred_sar[true, ], pch=1, cex=3,
         lwd=2)
  if(i == 1)
    legend('bottomright', 
           c('observed', 'iterative', 'iterative observed SAD', 'noniterative',
             'noniterative observed SAD'), pch=c(1, rep(NA, 4)), col=c(1, col),
           cex=2.5, bty='n', lwd=c(2, rep(4, 4)), lty=c(NA, lty))
}
dev.off()


## plot one-to-one curves-------------------------------------------------------

xlim = range(pred_sar[ , 1:4], na.rm=T)
ylim = range(pred_sar$richness)

jpeg('./figs/sar_one_to_one_mete_analytical.jpeg',
     width = 480 * 2, height= 480 * 2, quality = 100)
par(mfrow=c(2,2))
for(i in 1:4) {
  plot(richness ~ pred_sar[ , i], data = pred_sar, log='xy',
       axes=F, frame.plot=F, xlab='', ylab='',
       xlim=2^(c(-3, 9)), ylim=2^(c(-3, 9)))
  ticks = 2^seq(-3, 9, 2)
  addAxis1(at = ticks, lab = as.character(ticks))
  addAxis2(at = ticks, lab = as.character(ticks))
  lines(2^(c(-3, 9)), 2^(c(-3, 9)), lwd=2)
  r2 = get_R2(log(pred_sar$richness),
              log(pred_sar[ ,i]), na.rm=T)
  text(x=.2, y= 300, expression(italic(R^2) == ''), cex=3)
  text(x= 1, y = 275, round(r2, 3), cex=3)
  mtext(side=3, paste(letters[i], ')', sep=''), adj=0, cex=3)
}
dev.off()