## Author: Dan McGlinn
## Purpose: to create site x species matrices for the crosstimbers dataset
## which was surveyed in 1998 on a 64 m^2 plot.
## Metadata: CrosstimbersMasterDataSheet1-20-12.xlsx

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/cross1998_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$sp)))
dat$spnum = match(dat$sp, uniSpeciesNames)
S = max(dat$spnum) 

## for square quadrats the lengths are in meters defined here
i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0, 200, 0, 200) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$x, dat$y), n_quadrats, domain)

write.csv(comms, file='./data/cross_comms.csv', row.names=FALSE)