## METE analytical SAR predictions vs Empirical SARs

setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

fileNames = dir('./sar')

empirFiles = grep('empir_sar.csv', fileNames)
empir = vector('list', length(empirFiles))
names(empir) = sub('_empir_sar.csv', '', fileNames[empirFiles])
for (i in seq_along(empirFiles))
  empir[[i]] = read.csv(paste('./sar/', fileNames[empirFiles[i]], sep=''))

meteFiles = grep('mete_sar', fileNames)
mete = vector('list', length(meteFiles))
names(mete) = sub('_mete_sar.txt', '', fileNames[meteFiles])
for( i in seq_along(meteFiles)) {
  mete[[i]] = read.csv(paste('./sar/', fileNames[meteFiles[i]], sep=''))
  mete[[i]]$area = mete[[i]]$area * empir[[i]]$area[1]
}

## average cocoli and sherman plots
empir$sherman1 = (empir$sherman1 + empir$sherman2) / 2
empir$cocoli1 = (empir$cocoli1 + empir$cocoli2) / 2
mete$sherman1 = (mete$sherman1 + mete$sherman2) / 2
mete$cocoli1 = (mete$cocoli1 + mete$cocoli2) / 2
empir = empir[-match(c('sherman2','sherman3','cocoli2'), names(empir))]
mete = mete[-match(c('sherman2','sherman3','cocoli2'), names(mete))]

pdf('./figs/mete_&_empir_sar.pdf', width=7 * 2, height=7 * 2)
  par(mfrow=c(4,4))
  ## arithmetic
  for (i in seq_along(mete)) {
    plot(sr ~ area, data = mete[[i]], 
         ylim=range(c(mete[[i]]$sr, empir[[i]]$richness)),
         xlim=range(c(mete[[i]]$area, empir[[i]]$area)),
         type='o', main=names(mete)[i], xlab='Area (m2)')
    lines(richness ~ area, data = empir[[i]], pch=19, type='o')
    if(i == 1)
      legend('bottomright',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1)
  }
  ## log log
  for (i in seq_along(mete)) {
    plot(sr ~ area, data = mete[[i]], log='xy',
         ylim=range(c(mete[[i]]$sr, empir[[i]]$richness)),
         xlim=range(c(mete[[i]]$area, empir[[i]]$area)),
         type='o', main=names(mete)[i], xlab='Area (m2)')
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
  lines(sr ~ area, data = mete[[i]], pch=1, type='o', lwd=3, cex=2)
  lines(richness ~ area, data = empir[[i]], pch=19, type='o', lwd=2, cex=2)
}
par(mfrow=c(1,1))
plot(1:10,1:10,type='n',axes=F,xlab='',ylab='')
legend('center',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1,lwd=6,cex=3)

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
par(mfrow=c(1,1))
plot(predZ ~ logNS, data=slopes, type='n', ylim=c(0,1), axes=F, frame.plot=F,
     xlab='', ylab='')
axis(side=1,lwd=4,cex.axis=2)
axis(side=2,lwd=4,cex.axis=2)
lines(Zsmooth, lwd=4, lty=2)
comms = unique(slopes$comm)
col = colorRampPalette(c("blue", "red", "green", "orange"))(length(comms))
for (i in seq_along(comms))
  points(obsZ ~ logNS, data=slopes, subset=comm == comms[i], col=col[i], pch=19,
         cex=1.25)

par(mar=c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab='',ylab='')
legend('center',comms, pch=19, cex=2, col=col, bty='n')

    

