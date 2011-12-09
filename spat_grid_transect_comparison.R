##spatial summaries

pdf('grid_vs_transect.pdf',width=14,height=14)

setwd('/home/danmcglinn/maxent/trunk')
##varWithin

load('./varWithin/varWithin_S100_N10000_C200_B11_transect_binary.Rdata')
vTranExp = sapply(varWithin,function(x)x$varWithinObs$vario$exp.var)
distsTran = varWithin[[1]]$varWithinObs$vario$Dist
vTranExpAvg = apply(vTranExp,1,mean)
vTranExpQt = apply(vTranExp,1,function(x)quantile(x,c(.025,.975)))
##
load('./varWithin/varWithin_S100_N10000_C200_B11_grid_binary.Rdata')
vGridExp = sapply(varWithin,function(x)x$varWithinObs$vario$exp.var)
distsGrid = varWithin[[1]]$varWithinObs$vario$Dist
##scale 2D distances
distsGrid = distsGrid #* (max(distsTran)/max(distsGrid))
vGridExpAvg = apply(vGridExp,1,mean)
vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))

par(mfrow=c(2,2))
plot(distsTran,vTranExpAvg,type='n',ylim=range(vTranExpAvg,vTranExpQt),
     ylab='Variance Within-species',xlab='Spatial Lag',xlim=c(0,max(distsGrid)))
polygon(c(distsTran,distsTran[length(distsTran):1]),
        c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1]),col='pink',
        border = NA)
polygon(c(distsGrid,distsGrid[length(distsGrid):1]),
        c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1]),col='lightblue',
        border = NA)
lines(distsTran,vTranExpAvg,type='l',col='red',pch=19,lwd=2)
lines(distsGrid,vGridExpAvg,type='l',col='blue',pch=19,lwd=2)
legend("bottomright",c("Transect MaxEnt Avg","Transect MaxEnt CI",
       "Grid MaxEnt Avg","Grid MaxEnt CI"),cex=2,col=c('red','pink','blue',
       'lightblue'),lty=1,lwd=c(3,7,3,7),bty='n')
##
plot(log10(distsTran),log10(vTranExpAvg),type='n',ylim=log10(range(vTranExpAvg,
     vTranExpQt)),ylab='log10 Variance Within-species',xlab='log10 Spatial Lag',
     xlim=log10(range(distsGrid)))
polygon(log10(c(distsTran,distsTran[length(distsTran):1])),
        log10(c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1])),col='pink',
        border = NA)
polygon(log10(c(distsGrid,distsGrid[length(distsGrid):1])),
        log10(c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1])),
        col='lightblue',border = NA)
lines(log10(distsTran),log10(vTranExpAvg),type='l',col='red',pch=19,lwd=2)
lines(log10(distsGrid),log10(vGridExpAvg),type='l',col='blue',pch=19,lwd=2)

##sorensen
load('./sorensen/sorensen_S100_N10000_C200_B11_transect_binary.Rdata')
vTranExp = 1 - sapply(sorensen,function(x)x$sorensenObs$vario$exp.var)
distsTran = sorensen[[1]]$sorensenObs$vario$Dist
vTranExpAvg = apply(vTranExp,1,mean)
vTranExpQt = apply(vTranExp,1,function(x)quantile(x,c(.025,.975)))
##
load('./sorensen/sorensen_S100_N10000_C200_B11_grid_binary.Rdata')
vGridExp = 1 - sapply(sorensen,function(x)x$sorensenObs$vario$exp.var)
distsGrid = sorensen[[1]]$sorensenObs$vario$Dist
##scale 2D distances
distsGrid = distsGrid #* (max(distsTran)/max(distsGrid))
vGridExpAvg = apply(vGridExp,1,mean)
vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))

#par(mfrow=c(1,2))
plot(distsTran,vTranExpAvg,type='n',ylim=range(vTranExpAvg,vTranExpQt),
     ylab='Sorensen Dissimilarity',xlab='Spatial Lag',xlim=c(0,max(distsGrid)))
polygon(c(distsTran,distsTran[length(distsTran):1]),
        c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1]),col='pink',
        border = NA)
polygon(c(distsGrid,distsGrid[length(distsGrid):1]),
        c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1]),col='lightblue',
        border = NA)
lines(distsTran,vTranExpAvg,type='l',col='red',pch=19,lwd=2)
lines(distsGrid,vGridExpAvg,type='l',col='blue',pch=19,lwd=2)
legend("topright",c("Transect MaxEnt Avg","Transect MaxEnt CI",
       "Grid MaxEnt Avg","Grid MaxEnt CI"),cex=2,col=c('red','pink','blue',
       'lightblue'),lty=1,lwd=c(3,7,3,7),bty='n')
##
plot(log10(distsTran),log10(vTranExpAvg),type='n',ylim=log10(range(vTranExpAvg,
     vTranExpQt)),ylab='log10 Sorensen Dissimilarity',xlab='log10 Spatial Lag',
     xlim=log10(range(distsGrid)))
polygon(log10(c(distsTran,distsTran[length(distsTran):1])),
        log10(c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1])),col='pink',
        border = NA)
polygon(log10(c(distsGrid,distsGrid[length(distsGrid):1])),
        log10(c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1])),
        col='lightblue',border = NA)
lines(log10(distsTran),log10(vTranExpAvg),type='l',col='red',pch=19,lwd=2)
lines(log10(distsGrid),log10(vGridExpAvg),type='l',col='blue',pch=19,lwd=2)



### 2-D distances stretched to match linear distancs
load('./varWithin/varWithin_S100_N10000_C200_B11_transect_binary.Rdata')
vTranExp = sapply(varWithin,function(x)x$varWithinObs$vario$exp.var)
distsTran = varWithin[[1]]$varWithinObs$vario$Dist
vTranExpAvg = apply(vTranExp,1,mean)
vTranExpQt = apply(vTranExp,1,function(x)quantile(x,c(.025,.975)))
##
load('./varWithin/varWithin_S100_N10000_C200_B11_grid_binary.Rdata')
vGridExp = sapply(varWithin,function(x)x$varWithinObs$vario$exp.var)
distsGrid = varWithin[[1]]$varWithinObs$vario$Dist
##scale 2D distances
distsGrid = distsGrid * (max(distsTran)/max(distsGrid))
vGridExpAvg = apply(vGridExp,1,mean)
vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))

par(mfrow=c(2,2))
plot(distsTran,vTranExpAvg,type='n',ylim=range(vTranExpAvg,vTranExpQt),
     ylab='Variance Within-species',xlab='Spatial Lag',xlim=c(0,max(distsGrid)))
polygon(c(distsTran,distsTran[length(distsTran):1]),
        c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1]),col='pink',
        border = NA)
polygon(c(distsGrid,distsGrid[length(distsGrid):1]),
        c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1]),col='lightblue',
        border = NA)
lines(distsTran,vTranExpAvg,type='l',col='red',pch=19,lwd=2)
lines(distsGrid,vGridExpAvg,type='l',col='blue',pch=19,lwd=2)
legend("bottomright",c("Transect MaxEnt Avg","Transect MaxEnt CI",
       "Grid MaxEnt Avg","Grid MaxEnt CI"),cex=2,col=c('red','pink','blue',
       'lightblue'),lty=1,lwd=c(3,7,3,7),bty='n')
##
plot(log10(distsTran),log10(vTranExpAvg),type='n',ylim=log10(range(vTranExpAvg,
     vTranExpQt)),ylab='log10 Variance Within-species',xlab='log10 Spatial Lag',
     xlim=log10(range(distsGrid)))
polygon(log10(c(distsTran,distsTran[length(distsTran):1])),
        log10(c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1])),col='pink',
        border = NA)
polygon(log10(c(distsGrid,distsGrid[length(distsGrid):1])),
        log10(c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1])),
        col='lightblue',border = NA)
lines(log10(distsTran),log10(vTranExpAvg),type='l',col='red',pch=19,lwd=2)
lines(log10(distsGrid),log10(vGridExpAvg),type='l',col='blue',pch=19,lwd=2)

##sorensen
load('./sorensen/sorensen_S100_N10000_C200_B11_transect_binary.Rdata')
vTranExp = 1 - sapply(sorensen,function(x)x$sorensenObs$vario$exp.var)
distsTran = sorensen[[1]]$sorensenObs$vario$Dist
vTranExpAvg = apply(vTranExp,1,mean)
vTranExpQt = apply(vTranExp,1,function(x)quantile(x,c(.025,.975)))
##
load('./sorensen/sorensen_S100_N10000_C200_B11_grid_binary.Rdata')
vGridExp = 1 - sapply(sorensen,function(x)x$sorensenObs$vario$exp.var)
distsGrid = sorensen[[1]]$sorensenObs$vario$Dist
##scale 2D distances
distsGrid = distsGrid * (max(distsTran)/max(distsGrid))
vGridExpAvg = apply(vGridExp,1,mean)
vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))

#par(mfrow=c(1,2))
plot(distsTran,vTranExpAvg,type='n',ylim=range(vTranExpAvg,vTranExpQt),
     ylab='Sorensen Dissimilarity',xlab='Spatial Lag',xlim=c(0,max(distsGrid)))
polygon(c(distsTran,distsTran[length(distsTran):1]),
        c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1]),col='pink',
        border = NA)
polygon(c(distsGrid,distsGrid[length(distsGrid):1]),
        c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1]),col='lightblue',
        border = NA)
lines(distsTran,vTranExpAvg,type='l',col='red',pch=19,lwd=2)
lines(distsGrid,vGridExpAvg,type='l',col='blue',pch=19,lwd=2)
legend("topright",c("Transect MaxEnt Avg","Transect MaxEnt CI",
       "Grid MaxEnt Avg","Grid MaxEnt CI"),cex=2,col=c('red','pink','blue',
       'lightblue'),lty=1,lwd=c(3,7,3,7),bty='n')
##
plot(log10(distsTran),log10(vTranExpAvg),type='n',ylim=log10(range(vTranExpAvg,
     vTranExpQt)),ylab='log10 Sorensen Dissimilarity',xlab='log10 Spatial Lag',
     xlim=log10(range(distsGrid)))
polygon(log10(c(distsTran,distsTran[length(distsTran):1])),
        log10(c(vTranExpQt[1,],vTranExpQt[2,length(distsTran):1])),col='pink',
        border = NA)
polygon(log10(c(distsGrid,distsGrid[length(distsGrid):1])),
        log10(c(vGridExpQt[1,],vGridExpQt[2,length(distsGrid):1])),
        col='lightblue',border = NA)
lines(log10(distsTran),log10(vTranExpAvg),type='l',col='red',pch=19,lwd=2)
lines(log10(distsGrid),log10(vGridExpAvg),type='l',col='blue',pch=19,lwd=2)


dev.off()
