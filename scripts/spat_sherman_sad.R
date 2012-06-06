setwd('/home/danmcglinn/maxent/spat')

## read in comms file for the 1999 census (i.e., census 3) of sherman
comms = read.csv('./data/sherman_comms.csv',header=TRUE)
grains = unique(comms$grain)

## there is are two 200 x 100 m quadrats and one 140 x 140 m quad
grain1 = grains[grep('_1',grains)][1]
grain2 = grains[grep('_2',grains)][1]
grain3 = grains[grep('_3',grains)][1]

## calculate empirical SAD/RAD for each
sadAbs1 = apply(comms[comms$grain == grain1,-(1:3)],2,sum)
sadAbs2 = apply(comms[comms$grain == grain2,-(1:3)],2,sum)
sadAbs3 = apply(comms[comms$grain == grain3,-(1:3)],2,sum)

## drop species in each with abundances of zero
sadAbs1 = sadAbs1[sadAbs1>0]
sadAbs2 = sadAbs2[sadAbs2>0]
sadAbs3 = sadAbs3[sadAbs3>0]

write.table(matrix(sort(sadAbs1,dec=TRUE),nrow=1),
            file='./data/sherman1_sad.csv',sep=',',row.names=FALSE,
            col.names=FALSE)
write.table(matrix(sort(sadAbs2,dec=TRUE),nrow=1),
            file='./data/sherman2_sad.csv',sep=',',row.names=FALSE,
            col.names=FALSE)
write.table(matrix(sort(sadAbs3,dec=TRUE),nrow=1),
            file='./data/sherman3_sad.csv',sep=',',row.names=FALSE,
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

