##spatial summaries

pdf('directional_variograms.pdf',width=14,height=14)

setwd('/home/danmcglinn/maxent/trunk')
##varWithin

load('./varWithin/varWithin_S100_N10000_C200_B11_grid_0deg_binary.Rdata')
vExp = array(NA,c(2,length(varWithin[[1]]$varWithinObs$vario$exp.var),200))
vExp[1,,] = sapply(varWithin,function(x)x$varWithinObs$vario$exp.var)
load('./varWithin/varWithin_S100_N10000_C200_B11_grid_90deg_binary.Rdata')
vExp[2,,] = sapply(varWithin,function(x)x$varWithinObs$vario$exp.var)
dists = varWithin[[1]]$varWithinObs$vario$Dist
vExpAvg = apply(vExp,c(1,2),mean)
vExpQt = apply(vExp,c(1,2),function(x)quantile(x,c(.025,.975)))

par(mfrow=c(2,2))
plot(dists,vExpAvg[1,],type='n',ylim=range(vExpAvg,vExpQt),
     ylab='Variance Within-species',xlab='Spatial Lag')
for(i in 1:2){
#  polygon(c(dists,dists[length(dists):1]),
#          c(vExpQt[1,i,],vExpQt[2,i,length(dists):1]),col=i,
#          border = NA)
  lines(dists,vExpAvg[i,],type='o',col=i+2,pch=19,lwd=2)  
}
legend("bottomright",c('0 degrees','90 degrees'),cex=2,col=3:4,
       lwd=3,bty='n')
##
plot(log10(dists),log10(vExpAvg[1,]),type='n',ylim=log10(range(vExpAvg,
     vExpQt)),ylab='log10 Variance Within-species',xlab='log10 Spatial Lag')
for(i in 1:2){
#  polygon(log10(c(dists,dists[length(dists):1])),
#          log10(c(vExpQt[1,i,],vExpQt[2,i,length(dists):1])),col='pink',
#          border = NA)
  lines(log10(dists),log10(vExpAvg[i,]),type='o',col=i+2,pch=19,lwd=2)  
}

##
load('./sorensen/sorensen_S100_N10000_C200_B11_grid_0deg_binary.Rdata')
vExp = array(NA,c(2,length(sorensen[[1]]$sorensenObs$vario$exp.var),200))
vExp[1,,] = 1 - sapply(sorensen,function(x)x$sorensenObs$vario$exp.var)
load('./sorensen/sorensen_S100_N10000_C200_B11_grid_90deg_binary.Rdata')
vExp[2,,] = 1 - sapply(sorensen,function(x)x$sorensenObs$vario$exp.var)

dists = sorensen[[1]]$sorensenObs$vario$Dist
vExpAvg = apply(vExp,c(1,2),mean)
vExpQt = apply(vExp,c(1,2),function(x)quantile(x,c(.025,.975)))

plot(dists,vExpAvg[1,],type='n',ylim=range(vExpAvg,vExpQt),
     ylab='Sorensen Similarity',xlab='Spatial Lag')
for(i in 1:2){
#  polygon(c(dists,dists[length(dists):1]),
#          c(vExpQt[1,i,],vExpQt[2,i,length(dists):1]),col=i,
#          border = NA)
  lines(dists,vExpAvg[i,],type='o',col=i+2,pch=19,lwd=2)  
}
legend("topright",c('0 degrees','90 degrees'),cex=2,col=3:4,
       lwd=3,bty='n')
##
plot(log10(dists),log10(vExpAvg[1,]),type='n',ylim=log10(range(vExpAvg,
     vExpQt)),ylab='log10 Sorensen Similarity',xlab='log10 Spatial Lag')
for(i in 1:2){
#  polygon(log10(c(dists,dists[length(dists):1])),
#          log10(c(vExpQt[1,i,],vExpQt[2,i,length(dists):1])),col='pink',
#          border = NA)
  lines(log10(dists),log10(vExpAvg[i,]),type='o',col=i+2,pch=19,lwd=2)  
}

dev.off()
