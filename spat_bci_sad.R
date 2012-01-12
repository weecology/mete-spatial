setwd('/home/danmcglinn/maxent/spat')

## read in comms file for the 2010 census (i.e., census7) of BCI
comms = read.csv('./data/bci_comms.csv',header=TRUE)
grains = unique(comms[,1])

## calculate empirical SAD/RAD
sadAbs = apply(comms[comms[,1] == grains[1],-(1:3)],2,sum)
write.table(matrix(sort(sadAbs,dec=TRUE),nrow=1),
            file='./data/bci_sad.csv',sep=',',row.names=FALSE,
            col.names=FALSE)
sadAvg = sapply(grains,function(x) apply(comms[comms[,1]==x,-(1:3)],2,mean))

## graphically explore SAD patterns
S= ncol(comms)-3
par(mfrow=c(1,2))
plot(1:S,sort(sadAvg[,1],dec=TRUE),ylim=range(sadAvg),log='y',type='n')
for(i in seq_along(grains))
  lines(1:S,sort(sadAvg[,i]+1,dec=TRUE),type='o',col=i)
plot(1:S,sort(sadAvg[,1],dec=TRUE),ylim=range(sadAvg),log='xy',type='n')
for(i in seq_along(grains))
  lines(1:S,sort(sadAvg[,i]+1,dec=TRUE),type='o',col=i)


