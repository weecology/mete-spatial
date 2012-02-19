## METE analytical SAR predictions vs Empirical SARs

setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

AminExact = c(1e3/128 * 5e2/64,
              2e2/128 * 1e2/64,
              2e2/128 * 1e2/64,
              2e2/64 * 2e2/64,
              1,
              2e2/128 * 1e2/64,
              2e2/128 * 1e2/64,
              (1.4e2/64)^2)

fileNames = dir('./sar')
meteFiles = grep('mete_sar',fileNames)
mete = vector('list',length(meteFiles))
for(i in seq_along(meteFiles))
  mete[[i]] = read.csv(paste('./sar/',fileNames[meteFiles[i]],sep=''))
names(mete) = sub('_mete_sar.txt','',fileNames[meteFiles])

##for cocoli1, cocoli2, sherman1, sherman2 we need to adjust the sars
##b/c the mete prediction fails at scales where there is less than 1 individual
meteNames = names(mete)
for(i in seq_along(mete)){
  if(any(meteNames[i] == c('cocoli1', 'cocoli2', 'cross','sherman1', 'sherman2'))){
    mete[[i]]$area = mete[[i]]$area * 2
    mete[[i]] = rbind(c(1,NA),mete[[i]])
  }
  mete[[i]]$area = mete[[i]]$area * AminExact[i]
}

empirFiles = grep('empir',fileNames)
empir = vector('list',length(empirFiles))
for(i in seq_along(empirFiles)){
  empir[[i]] = read.csv(paste('./sar/',fileNames[empirFiles[i]],sep=''))
  empir[[i]] = empir[[i]][,-3]
}  
names(empir) = sub('_empir_sar.csv','',fileNames[empirFiles])



pdf('METE_&_empir_SAR.pdf',width=14,height=7)
par(mfrow=c(2,4))
for(i in seq_along(mete)){
  plot(mete[[i]],ylim=range(c(mete[[i]][,2],empir[[i]][,2]),na.rm=TRUE),
       xlim=range(c(mete[[i]][,1],empir[[i]][,1])),type='o',
       main=names(mete)[i],xlab='Area (m2)')
  lines(empir[[i]],pch=19,type='o')
  if(i == 1)
    legend('bottomright',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1)
}
  
par(mfrow=c(2,4))
for(i in seq_along(mete)){
  plot(mete[[i]],ylim=range(c(mete[[i]][,2],empir[[i]][,2]),na.rm=TRUE),
       xlim=range(c(mete[[i]][,1],empir[[i]][,1])),type='o',log='xy',
       main=names(mete)[i],xlab='Area (m2)')
  lines(empir[[i]],pch=19,type='o')
  if(i == 1)
    legend('bottomright',c('Empirical','METE'),pch=c(19,1),bty='n',lty=1)
}
dev.off()    
  
## for powerpoint
par(mfrow=c(2,4))
for(i in seq_along(mete)){
  plot(mete[[i]],ylim=range(c(mete[[i]][,2],empir[[i]][,2]),na.rm=TRUE),
       xlim=range(c(mete[[i]][,1],empir[[i]][,1])),type='n',
       main=names(mete)[i],ylab='',xlab='',axes=F)
  axis(side=1,lwd=4,cex.axis=2)
  axis(side=2,lwd=4,cex.axis=2)
  lines(mete[[i]],pch=1,type='o',lwd=3,cex=2)
  lines(empir[[i]],pch=19,type='o',lwd=2,cex=2)
}
  
par(mfrow=c(2,4))
for(i in seq_along(mete)){
  plot(mete[[i]],ylim=range(c(mete[[i]][,2],empir[[i]][,2]),na.rm=TRUE),
       xlim=range(c(mete[[i]][,1],empir[[i]][,1])),type='n',
       main=names(mete)[i],ylab='',xlab='',axes=F,log='xy')
  axis(side=1,lwd=4,cex.axis=2)
  axis(side=2,lwd=4,cex.axis=2)
  lines(mete[[i]],pch=1,type='o',lwd=3,cex=2)
  lines(empir[[i]],pch=19,type='o',lwd=3,cex=2)
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
par(mfrow=c(1,1))
plot(predZ ~ logNS,data=slopes,type='n',ylim=c(0,1),axes=F,frame.plot=F,xlab='',ylab='')
axis(side=1,lwd=4,cex.axis=2)
axis(side=2,lwd=4,cex.axis=2)
points(predZ ~ logNS,data=slopes,col='black',pch=1,cex=2)
cls = c('red',rep('blue',2),'dodgerblue','green3',rep('purple3',3))
datNames = c('bci','cocoli1','cocoli2','cross','serp','sherman1','sherman2','sherman3')
for(i in 1:length(datNames))
  points(obsZ ~ logNS,data=slopes,subset=comm==datNames[i],col=cls[i],pch=19,cex=1.25)

par(mar=c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab='',ylab='')
legend('center',c('METE','BCI','Cocoli','Crosstimbers','Serpentine','Sherman'),
       pch=c(1,rep(19,5)),cex=2,col=c('black',unique(cls)),bty='n')

    

