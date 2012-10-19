## METE analytical SAR predictions vs Empirical SARs

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

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

## load expected sars under random placement
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

## plot an example SAR figure for presentation
par(mfrow=c(1,3))
site = 'sherman1'
purple = rgb(112, 48, 160, maxColorValue=160)
lightblue = "#1AB2FF"

i = match(site, names(mete))
    plot(sr ~ area, data=mete[[i]], log='xy',
         xlim=c(1e1, 1e6), ylim=c(1e1, 1e3),
         type='n', xlab='', ylab='',frame.plot=F, axes=F,
         main = site)
    axis(side=1, at=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6))
    axis(side=2, at=c(1e1, 1e2, 1e3))
    ## mete CI
#    dat = meteavg[[match(names(mete)[i], names(meteavg))]]
#    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
#    lines(sr.avg ~ grains, data=dat, col=lightblue, lwd=4, lty=2)
    ## RP CI
    dat = srExp[[match(names(mete)[i], names(srExp))]]
#    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='pink')
    lines(srCol~ grains, data=dat, type='p', pch=15, col='blue', lwd=4,
          lty=2, cex=2)
    ## analytical mete    
    lines(sr ~ area, data=mete[[i]], type='p', pch=15, col='green3',
          lwd=4, cex=2)
    ## data
    lines(richness ~ area, data = empir[[i]], pch=0, type='p',
          lwd=2, cex=3)


pdf('./figs/mete_&_empir_sar.pdf', width=7 * 2, height=7 * 2)
  par(mfrow=c(4,4))
  ## arithmetic
  for (i in seq_along(mete)) {
    plot(sr ~ area, data=mete[[i]], 
         ylim=range(c(mete[[i]]$sr, empir[[i]]$richness)),
         xlim=range(c(mete[[i]]$area, empir[[i]]$area)),
         type='n', main=names(mete)[i], xlab='Area (m2)')
    ## mete CI
    if (names(mete)[i] != 'cross') {
      dat = meteavg[[match(names(mete)[i], names(meteavg))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(sr.avg ~ grains, data=dat, col='red', type='o')
    }  
    ## RP CI
    dat = srExp[[match(names(mete)[i], names(srExp))]]
    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='pink')
    lines(srCol~ grains, data=dat, col='red', type='o')
    ## analytical mete    
    lines(sr ~ area, data=mete[[i]], type='o')
    ## data
    lines(richness ~ area, data = empir[[i]], pch=19, type='o')
    if(i == 1)
      legend('bottomright',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1)
  }
  ## log log
  for (i in seq_along(mete)) {
    plot(sr ~ area, data=mete[[i]], log='xy',
         ylim=range(c(mete[[i]]$sr, empir[[i]]$richness)),
         xlim=range(c(mete[[i]]$area, empir[[i]]$area)),
         type='n', main=names(mete)[i], xlab='Area (m2)')
    if (names(mete)[i] != 'cross') {
      dat = meteavg[[match(names(mete)[i], names(meteavg))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(sr.avg ~ grains, data=dat, col='red', type='o')
    }
    ## RP CI
    dat = srExp[[match(names(mete)[i], names(srExp))]]
    addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='pink')
    lines(srCol~ grains, data=dat, col='red', type='o')
    ## analytical mete    
    lines(sr ~ area, data=mete[[i]], type='o')
    ## data
    lines(richness ~ area, data = empir[[i]], pch=19, type='o')
    if(i == 1)
      legend('bottomright',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1)
  }
dev.off()    
  
## for powerpoint
## arith
par(mfrow=c(4,4))
for (i in seq_along(mete)) {
  plot(sr ~ area, data = mete[[i]], 
       ylim=range(c(mete[[i]]$sr, empir[[i]]$richness)),
       xlim=range(c(mete[[i]]$area, empir[[i]]$area)),
       type='n', main=names(mete)[i], ylab='',xlab='',axes=F)  
  axis(side=1, lwd=4, cex.axis=2)
  axis(side=2, lwd=4, cex.axis=2)
  if (names(mete)[i] != 'cross') {
      dat = meteavg[[match(names(mete)[i], names(meteavg))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(sr.avg ~ grains, data=dat, col='red', type='o')
  }      
  lines(sr ~ area, data = mete[[i]], pch=1, type='o', lwd=3, cex=2)
  lines(richness ~ area, data = empir[[i]], pch=19, type='o', lwd=2, cex=2)
}
## loglog
par(mfrow=c(4,4))
for (i in seq_along(mete)) {
  plot(sr ~ area, data = mete[[i]], log='xy',
       ylim=range(c(mete[[i]]$sr, empir[[i]]$richness)),
       xlim=range(c(mete[[i]]$area, empir[[i]]$area)),
       type='n', main=names(mete)[i], ylab='',xlab='',axes=F)  
  axis(side=1, lwd=4, cex.axis=2)
  axis(side=2, lwd=4, cex.axis=2)
  if (names(mete)[i] != 'cross') {
      dat = meteavg[[match(names(mete)[i], names(meteavg))]]
      addCI('grains', 'sr.lo', 'sr.hi', data='dat', col='grey')
      lines(sr.avg ~ grains, data=dat, col='red', type='o')
  }        
  lines(sr ~ area, data = mete[[i]], pch=1, type='o', lwd=3, cex=2)
  lines(richness ~ area, data = empir[[i]], pch=19, type='o', lwd=2, cex=2)
}
par(mfrow=c(1,1))
plot(1:10,1:10,type='n',axes=F,xlab='',ylab='')
legend('center',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1,lwd=6,cex=3)

## compare residuals between METE and RP
## create flat file for easy access

for (i in seq_along(empir)) {
  site = names(empir)[i]
  if (site == 'cross')
    next
  area = empir[[i]]$area
  obs.sr = empir[[i]]$richness
  index = match(site, names(meteavg))
  mete.sr = meteavg[[index]]$sr.avg
  mete.res = obs.sr - mete.sr
  index = match(site, names(srExp))
  rp.sr = srExp[[index]]$srCol
  rp.res = obs.sr - rp.sr
  if (exists('sr_res'))
    sr_res = rbind(sr_res,
                   data.frame(site, area, obs.sr, mete.sr, rp.sr, mete.res, rp.res))
  else
    sr_res = data.frame(site, area, obs.sr, mete.sr, rp.sr, mete.res, rp.res)
}

plot(mete.sr ~ obs.sr, data=sr_res)
points(rp.sr ~ obs.sr, data=sr_res, pch=19)
abline(a=0, b=1)

plot(mete.res ~ rp.res, data=sr_res)
abline(a=0, b=1)

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

##----------------------------------------------------------------------------
## plot universal curve

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
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
hab = c('tropical', 'oak-hickory', 'pine', 'oak savanna', 'mixed evergreen',
        'grassland')
col = c("forestgreen", "#1AB2FF", "medium purple", "#E61A33", "#B2FF8C",
        "#FF8000")
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

##------------------------------------------------------------------------------
## examine empirical SAR and averaged METE sar 
#bring in serp mete avg and quantile

serp = read.csv('./sar/serp_mete_avg.csv')
par(mfrow=c(1,1))
plot(mete$serp,type='n',log='xy',ylim=range(list(mete$serp$sr,serp[,2:4])),
     frame.plot=F,axes=F,xlab='',ylab='')
axis(side=1,lwd=4,cex.axis=2)
axis(side=2,lwd=4,cex.axis=2)
polygon(c(serp$area,rev(serp$area)),c(serp$sLow,rev(serp$sHigh)),border=NA,col='grey')
points(serp[,1:2],type='l',lwd=2)
points(mete$serp,cex=1.25)
points(empir$serp,pch=19,cex=1.25)
legend('bottomright',c('Empirical','METE analytical','METE simlulated'),
       bty='n',lty=c(NA,NA,1),lwd=c(NA,NA,4),pch=c(1,19,NA),cex=1.5)
       

