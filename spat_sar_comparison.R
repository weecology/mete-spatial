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
meteFiles = grep('mete',fileNames)
mete = vector('list',length(meteFiles))
for(i in seq_along(meteFiles))
  mete[[i]] = read.csv(paste('./sar/',fileNames[meteFiles[i]],sep=''))
names(mete) = sub('_mete_sar.txt','',fileNames[meteFiles])

##for cocoli1, cocoli2, sherman1, sherman2 we need to adjust the sars
##b/c the mete prediction fails at scales where there is less than 1 individual
toFix = c('cocoli1', 'cocoli2', 'sherman1', 'sherman2')
toFix = which(names(mete)%in%toFix)
for(i in seq_along(mete)){
  if(i %in% toFix){
    mete[[i]]$area = mete[[i]]$area * 2
    mete[[i]] = rbind(c(1,NA),mete[[i]])
  }
  mete[[i]]$area = mete[[i]]$area * AminExact[i]
}
toFix = 'cross'
toFix = which(names(mete)%in%toFix)
for(i in seq_along(toFix)){
  mete[[toFix[i]]]$area = mete[[toFix[i]]]$area * 4
  mete[[toFix[i]]] = rbind(c(1,NA),c(2,NA),mete[[toFix[i]]])
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

