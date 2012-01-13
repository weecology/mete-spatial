## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 7
## which was surveyed in 2010 for each spatial scale of interest
## Metadata: bci50ha.doc 

source('/home/danmcglinn/maxent/spat/spat_sim_vario_func.R')

## read in data from the 2010 census (i.e. census 7)

load('/home/danmcglinn/CTFSplots/BCI/bci_census7.Rdata')

uniSpeciesNames = as.character(sort(unique(dat$Latin)))
dat$spnum = match(dat$Latin,uniSpeciesNames)
S = max(dat$spnum) 

shortSide = 500
longSide = 1000

## for square quadrats the lengths are in meters defined here
nPixels = wid(c(14,12,10,8,6,4))
quadLen = shortSide/ nPixels
quadN = nPixels^2 * (longSide / shortSide)

## generate a site x species matrix for each spatial scale

comms = makeCommMat(dat$spnum,S,cbind(dat$gx,dat$gy),quadLen,quadN,c(0,1000,0,500))

write.csv(comms,file='/home/danmcglinn/maxent/spat/data/bci_comms.csv',
          row.names=FALSE)
save(comms,file='/home/danmcglinn/CTFSplots/BCI/bci_comms_census7.Rdata')
     