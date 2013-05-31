setwd('~/maxent/spat')

source('./scripts/spat_functions.R')

dir.create('./figs')

print('Generating SAR figures, ...')

## load data
sar_res = read.csv('./sar/sar_residuals.csv')

## drop all rows in which density of indivdiuals is less than 2
sar_res = sar_res[sar_res$indiv > 2, ]

indices = c(3:4, 6:7)

pred_sar = sar_res$richness - sar_res[ , indices]
pred_sar = pred_sar[ , c(1, 3, 2, 4)]
pred_sar = data.frame(pred_sar, richness = sar_res$richness, 
                      area = sar_res$area, site = sar_res$site,
                      hab = sar_res$hab)

## plot individual site relationships-------------------------------------------

site_names = as.character(as.matrix(read.table('./data/shrtnames.txt')))
site_titles = sub('1', '', site_names)

capwords = function(s, strict = FALSE) {
  cap = function(s) paste(toupper(substring(s, 1, 1)),
                          {s = substring(s, 2); if(strict) tolower(s) else s},
                          sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

site_titles = capwords(site_titles)

col = c('red', 'red', 'dodgerblue', 'dodgerblue')
lty = c(1, 3, 1, 3)

## log-log plots
jpeg('./figs/mete_&_empir_sar_log2.jpeg', width=480 * 4, height=480 * 4, quality=100)
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
    mtext(side=3, paste(site_titles[i], '-', unique(pred_sar$hab[true])),
          cex=2)
    mtext(side=3, paste('(', letters[i], ')', sep=''), adj=0, cex=2, font=2)
    for (j in 1:4) 
      lines(log2(pred_sar[true, j]) ~ log2(area), data=pred_sar[true, ], col=col[j],
            lwd=4, lty=lty[j])
    ## data
    points(log2(richness) ~ log2(area), data = pred_sar[true, ], pch=1, cex=3,
           lwd=2)
    if(i == 1)
      legend('bottomright', 
             c('observed', 'recursive', 'recursive, observed SAD', 'non-recursive',
               'non-recursive, observed SAD'), pch=c(1, rep(NA, 4)), col=c(1, col),
             cex=2.5, bty='n', lwd=c(2, rep(4, 4)), lty=c(NA, lty))
  }
dev.off()


## plot one-to-one curves-------------------------------------------------------

titles = c('recursive', 'recursive, observed SAD', 'non-recursive',
           'non-recursive, observed SAD')

jpeg('./figs/sar_one_to_one_mete_analytical.jpeg',
     width = 480 * 2, height= 480 * 2, quality = 100)
  par(mfrow=c(2,2))
  for(i in 1:4) {
    plot(richness ~ pred_sar[ , i], data = pred_sar, log='xy',
         axes=F, frame.plot=F, xlab='', ylab='',
         xlim=2^(c(-1, 9)), ylim=2^(c(-1, 9)))
    ticks = 2^seq(-1, 9, 2)
    addAxis(side=1, at = ticks, lab = as.character(ticks))
    addAxis(side=2, at = ticks, lab = as.character(ticks), padj=0)
    lines(2^(c(-1, 9)), 2^(c(-1, 9)), lwd=2)
    r2 = get_R2(log(pred_sar$richness), log(pred_sar[ , i]), na.rm=T)
    r2 = as.character(round(r2, 3))
    if (nchar(r2) < 5)
      r2 = paste(r2, rep(0, 5 - nchar(r2)), sep='')
    text(x=.9, y= 300, expression(italic(R^2) == ''), cex=2.5)
    text(x=2.75, y= 275, r2, cex=2.5)
    mtext(side=3, titles[i], cex=2)
    mtext(side=3, paste('(', letters[i], ')', sep=''), adj=0, cex=2, font=2)
  }
dev.off()

print('Generating SAR figures, complete!')

