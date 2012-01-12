## Purpose: to analyze Green et al's serpentine dataset
## data here: http://socrates.berkeley.edu/~hartelab/articles/serpentine_data.xls
## metadata here: http://socrates.berkeley.edu/~hartelab/MaxEnt.html
setwd('/home/danmcglinn/maxent')

library(vegan)
library(danspkg)
source('spat_sim_vario_func.R')

dat = read.csv('./data/serpentine_data.csv',header=T)
dim(dat)
# 259 x 24
head(dat)
tail(dat)
## drop last three garbage rows
dat = dat[-(nrow(dat):(nrow(dat)-2)),]
tail(dat)
## fix random NA in row 1
dat[1,] = ifelse(is.na(dat[1,]),0,dat[1,])
## calculate S and N
ncol(dat) # S = 24
sum(dat) # N = 37182

x = rep(1:16,each=16)
y = rep(1:16,times=16)
grains = c(1,4,16)
comms = quadAggregator(mat = dat, coord = cbind(x,y), grains = grains)

## calculate empirical Distance Decay patterns
metricsAbu = calcMetrics(comms,'all','abu',writeToFile=TRUE,fileSuffix='_serp')
metricsBinary = calcMetrics(comms,'all','binary',writeToFile=TRUE,fileSuffix='_serp')

