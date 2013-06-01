setwd('~/maxent/spat')
source('./scripts/spat_functions.R')

print('Generating SAR figures, ...')

## load data
sar_res = read.csv('./sar/sar_residuals.csv')

## drop all rows in which density of indivdiuals is less than 2
sar_res = sar_res[sar_res$indiv > 2, ]

indices = c(4:5, 8:9)

pred_sar = sar_res$richness - sar_res[ , indices]
pred_sar = pred_sar[ , c(1, 3, 2, 4)]
pred_sar = data.frame(pred_sar, richness = sar_res$richness, 
                      area = sar_res$area, site = sar_res$site,
                      hab = sar_res$hab)

## plot individual site relationships-------------------------------------------

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
#site_titles[1] = "BCI"
#site_titles[11] = "Bald Mtn."
#site_titles[14] = "UCSC"
#site_titles[15] = "Serpentine"
#site_titles[16] = "Cross Timbers"

col = c('red', 'red', 'dodgerblue', 'dodgerblue')
lty = c(1, 3, 1, 3)

## arthimetic plots
jpeg('./figs/mete_&_empir_sar.jpeg', width=480 * 4, height=480 * 4, quality=100)
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
    mtext(side=3, paste(site_titles[i], '-', unique(pred_sar$hab[true])),
          cex=2)
    mtext(side=3, paste('(', letters[i], ')', sep=''), adj=0, cex=2, font=2)
    for (j in 1:4) 
      lines(pred_sar[true, j] ~ area, data=pred_sar[true, ], col=col[j], lwd=4,
            lty=lty[j])
    ## data
    points(richness ~ area, data = pred_sar[true, ], pch=1, cex=3, lwd=2)
    if(i == 1)
      legend('bottomright', 
             c('observed', 'recursive', 'recursive, observed SAD', 'non-recursive',
               'non-recursive, observed SAD'), pch=c(1, rep(NA, 4)), col=c(1, col),
             cex=2.5, bty='n', lwd=c(2, rep(4, 4)), lty=c(NA, lty))
  }
dev.off()

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

## supplemental cocoli and sherman study site figure ---------------------------
mk_sup_figs = FALSE
if (mk_sup_figs) {
  cocoli = read.csv('./data/filtered_data/cocoli_census3_filtered.csv')
  sherman = read.csv('./data/filtered_data/sherman_census3_filtered.csv')
  ucsc = read.csv('./data/filtered_data/ferp_2007_filtered.csv')
  luquillo = read.csv('./data/filtered_data/luquillo_census4_filtered.csv')
  
  
  ## split plot into two 200 x 100 m rectangles
  cocoli1 = cocoli[cocoli$y>=100,]
  cocoli2 = cocoli[cocoli$y<100,]
  
  ## split plot into three quadrats (two 200 x 100 m rectangles and one 140 x 140m)
  sherman1 = sherman[sherman$x >= 140 & sherman$y >= 240, ]  ## one of the 200 x 100 m rectangles
  sherman2 = sherman[sherman$x >= 140 & sherman$y < 240, ] ## the other 200 x 100 m rectangle
  sherman3 = sherman[sherman$x < 140, ]  ## the 140 x 140 square
  
  cex = 2
  
  jpeg('./figs/subplot_maps.jpeg',
       width=480 * 2, height=480 * 2, quality=100)
  par(mfrow=c(2,2))
  plot(cocoli1$x,cocoli1$y,col='red', pch='.', ylim=range(sherman$y),
       xlim=range(sherman$x), xlab='X-coordinate', ylab='Y-coordinate',
       axes=F, frame.plot=T, cex=cex)
  points(cocoli2$x,cocoli2$y,col='blue',pch='.', cex=cex)
  lines(c(0, 0), c(0, 300), lwd=2)
  lines(c(0, 200), c(0, 0), lwd=2)
  lines(c(200, 200), c(0, 100), lwd=2)
  lines(c(0, 200), c(100, 100), lwd=2)
  lines(c(100, 100), c(100, 300), lwd=2)
  lines(c(0, 100), c(300, 300), lwd=2)
  axis(side=1, at = c(0, 100, 200))
  axis(side=2)
  mtext(side=3, 'Cocoli', cex=2.5)
  legend('topleft', c('Cocoli subplot 1', 'Cocoli subplot 2'),
         col=c('red','blue'), pch=19, bty='n', cex=cex)
  ##
  plot(sherman1$x,sherman1$y,col='red', pch='.', ylim=range(sherman$y),
       xlim=range(sherman$x), xlab='X-coordinate', ylab='Y-coordinate',
       axes=F, frame.plot=T, cex=cex)
  points(sherman2$x,sherman2$y,col='blue',pch='.', cex=cex)
  points(sherman3$x,sherman3$y,col='grey',pch='.', cex=cex)
  lines(c(140, 240), c(40, 40), lwd=2)
  lines(c(140, 140), c(40, 440), lwd=2)
  lines(c(140, 240), c(440, 440), lwd=2)
  lines(c(240, 240), c(40, 440), lwd=2)
  lines(c(140, 240), c(240, 240), lwd=2)
  axis(side=1, at = c(0, 100, 200))
  axis(side=2)
  mtext(side=3, 'Sherman', cex=2.5)
  legend('topleft', c('Sherman subplot 1', 'Sherman subplot 2'),
         col=c('red','blue'), pch=19, bty='n', cex=cex)
  ##
  trim = (320 - 250) / 2
  domain = c(trim, 320 - trim, 0, 500) # spatial domain in meters defined here
  true = luquillo$GX >= domain[1] & luquillo$GX < domain[2]
  plot(luquillo$GX[!true], luquillo$GY[!true], col='grey', pch='.', xlab='X-coordinate',
       ylab='Y-coordinate', xlim=range(luquillo$GX), ylim=range(luquillo$GY),
       axes=F, frame.plot=T, cex=cex)
  points(luquillo$GX[true], luquillo$GY[true], col='red', pch='.')
  lines(c(trim, 320 - trim), c(0, 0), lwd=2)
  lines(c(trim, trim), c(0, 500), lwd=2)
  lines(c(trim, 320 - trim), c(500, 500), lwd=2)
  lines(c(320 - trim, 320 - trim), c(0, 500), lwd=2)
  axis(side=1, c(0, 100, 200, 300))
  axis(side=2, c(0, 100, 200, 300, 400, 500))
  mtext(side=3, 'Luquillo', cex=2.5)
  ## 
  trim = (200 - 150) / 2
  domain = c(trim, 200 - trim, 0, 300) # spatial domain in meters defined here
  true = ucsc$gx >= domain[1] & ucsc$gx < domain[2]
  plot(ucsc$gx[!true], ucsc$gy[!true], col='grey', pch='.', xlab='X-coordinate',
       ylab='Y-coordinate', xlim=range(ucsc$gx), ylim=range(ucsc$gy),
       axes=F, frame.plot=T, cex=cex)
  points(ucsc$gx[true], ucsc$gy[true], col='red', pch='.')
  lines(c(trim, 200 - trim), c(0, 0), lwd=2)
  lines(c(trim, trim), c(0, 300), lwd=2)
  lines(c(trim, 200 - trim), c(300, 300), lwd=2)
  lines(c(200 - trim, 200 - trim), c(0, 300), lwd=2)
  axis(side=1, c(0, 100, 200))
  axis(side=2, c(0, 100, 200, 300))
  mtext(side=3, 'UCSC', cex=2.5)
  
  dev.off()
}

print('Generating SAR figures, complete!.')

