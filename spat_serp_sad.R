setwd('/home/danmcglinn/maxent/spat')

comms = read.csv('./data/serp_comms.csv')
grains = unique(comms$grain)

## calculate empirical SAD/RAD
sadAbs = apply(comms[comms[,1]==grains[1],-(1:3)],2,sum)
write.table(matrix(sort(sadAbs,dec=TRUE),nrow=1),
            file='./data/serpentine_sad.csv',sep=',',row.names=FALSE,
            col.names=FALSE)
sadAvg = sapply(grains,function(x) apply( comms[comms[,1]==x,-(1:3)],2,mean))
radAvg = sadAvg/apply(sadAvg,2,sum)

## graphically explore SAD/RAD patterns
S = ncol(comms)-3
par(mfrow=c(1,2))
plot(1:S,sort(sadAvg[,1],dec=TRUE),ylim=range(sadAvg),log='xy',type='n')
for(i in seq_along(grains))
  lines(1:S,sort(sadAvg[,i],dec=TRUE),type='o',col=i)
plot(1:S,sort(radAvg[,1],dec=TRUE),ylim=range(radAvg),log='xy',type='o')
for(i in seq_along(grains))
  lines(1:S,sort(radAvg[,i],dec=TRUE),type='o',col=i)
