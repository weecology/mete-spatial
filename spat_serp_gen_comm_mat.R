## Author: Dan McGlinn
## Purpose: to create site x species matrices for the serpentine dataset
## which was surveyed in 1998 on a 64 m^2 plot.
## Metadata: http://socrates.berkeley.edu/~hartelab/MaxEnt.html

setwd('/home/danmcglinn/maxent/spat')

source('spat_sim_vario_func.R')

dat = read.csv('./data/serpentine_data.csv',header=T)
dat = dat[-(nrow(dat):(nrow(dat)-2)),]
dat[1,] = ifelse(is.na(dat[1,]),0,dat[1,])
x = rep(1:16,each=16)
y = rep(1:16,times=16)
grains = c(1,4,16)
comms = quadAggregator(mat = dat, coord = cbind(x,y), grains = grains)

write.csv(comms,file='./data/serp_comms.csv',row.names=FALSE)