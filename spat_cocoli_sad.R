setwd('/home/danmcglinn/maxent/spat')

## read in comms file for the 1998 census (i.e., census 3) of cocoli
comms = read.csv('./data/cocoli_comms.csv',header=TRUE)
grains = unique(comms[,1])

## there is two 200 x 100 m quadrats

## calculate empirical SAD/RAD for each
sadAbs1 = apply(comms[comms[,1] == grains[1],-(1:3)],2,sum)
sadAbs2 = apply(comms[comms[,1] == grains[7],-(1:3)],2,sum)

## drop species in each with abundances of zero
sadAbs1 = sadAbs1[sadAbs1>0]
sadAbs2 = sadAbs1[sadAbs2>0]

write.table(matrix(sort(sadAbs1,dec=TRUE),nrow=1),
            file='./data/cocoli_sad_1.csv',sep=',',row.names=FALSE,
            col.names=FALSE)
write.table(matrix(sort(sadAbs2,dec=TRUE),nrow=1),
            file='./data/cocoli_sad_2.csv',sep=',',row.names=FALSE,
            col.names=FALSE)
            
sadAvg = sapply(grains,function(x) apply(comms[comms[,1]==x,-(1:3)],2,mean))

## graphically explore SAD patterns
S= ncol(comms)-3
par(mfrow=c(1,2))
plot(1:S,sort(sadAvg[,1]+1,dec=TRUE),ylim=range(sadAvg)+1,log='y',type='n')
for(i in seq_along(grains))
  lines(1:S,sort(sadAvg[,i]+1,dec=TRUE),type='o',col=i)
plot(1:S,sort(sadAvg[,1]+1,dec=TRUE),ylim=range(sadAvg)+1,log='xy',type='n')
for(i in seq_along(grains))
  lines(1:S,sort(sadAvg[,i]+1,dec=TRUE),type='o',col=i)

