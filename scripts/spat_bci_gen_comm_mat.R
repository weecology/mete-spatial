## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 7
## which was surveyed in 2010 for each spatial scale of interest
## Metadata: bci50ha.doc 

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

## read in data from the 2010 census (i.e. census 7)

load('~/datasets/CTFSplots/BCI/bci_census7.Rdata')

uniSpeciesNames = as.character(sort(unique(dat$Latin)))
dat$spnum = match(dat$Latin,uniSpeciesNames)
S = max(dat$spnum) 

i_bisections = c(13, 11, 9, 7, 5, 3)
n_quadrats = 2^i_bisections
domain = c(0,1000,0,500) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$gx, dat$gy), n_quadrats, domain)

write.csv(comms, file='./data/bci_comms.csv', row.names=FALSE)