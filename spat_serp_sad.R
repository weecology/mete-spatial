setwd('/home/danmcglinn/maxent')

source('spat_sim_vario_func.R')

dat = read.csv('./data/serpentine_data.csv',header=T)
dat = dat[-(nrow(dat):(nrow(dat)-2)),]
dat[1,] = ifelse(is.na(dat[1,]),0,dat[1,])
x = rep(1:16,each=16)
y = rep(1:16,times=16)
grains = c(1,4,16)
comms = quadAggregator(mat = dat, coord = cbind(x,y), grains = grains)

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
