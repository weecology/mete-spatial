## Author: Dan McGlinn
## Purpose: to create site x species matrices for the oosting tree dataset
## Metadata: http://esapubs.org/archive/ecol/E088/162/metadata.htm

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/oosting_trees_1990_filtered.csv')

S = length(unique(dat$CODE))
dat$spnum = as.integer(dat$CODE)

i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain = c(0,160,0,160)

comms = make_comm_matrix(dat$spnum, S, cbind(dat$NX, dat$NY), n_quadrats, domain)

## output results
write.csv(comms,file='./data/oosting_comms.csv',row.names=F)