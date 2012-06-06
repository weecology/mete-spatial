## Author: Dan McGlinn
## Purpose: to create site x species matrices for the crosstimbers dataset
## which was surveyed in 1998 on a 64 m^2 plot.
## Metadata: CrosstimbersMasterDataSheet1-20-12.xlsx

setwd('/home/danmcglinn/maxent/spat')

source('spat_sim_vario_func.R')

dat = read.csv('/home/danmcglinn/datasets/crosstimbers/cross1998_filtered.csv')
uniSpeciesNames = as.character(sort(unique(dat$sp)))
dat$spnum = match(dat$sp,uniSpeciesNames)
S = max(dat$spnum) 

shortSide = 200
longSide = 200

## for square quadrats the lengths are in meters defined here
nPixels = wid(c(14,12,10,8,6,4))
quadLen = shortSide/ nPixels
quadN = nPixels^2 * (longSide / shortSide)

## generate a site x species matrix for each spatial scale

comms = makeCommMat(dat$spnum,S,cbind(dat$x,dat$y),quadLen,quadN,c(0,200,0,200))

write.csv(comms,file='/home/danmcglinn/maxent/spat/data/cross_comms.csv',
          row.names=FALSE)