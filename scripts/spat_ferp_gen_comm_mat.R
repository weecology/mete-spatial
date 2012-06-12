
setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/ferp_2007_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$Latin)))
dat$spnum = match(dat$Latin, uniSpeciesNames)
S = max(dat$spnum) 

range(dat$gx) ## max 200
range(dat$gy) ## max 300
## change into a single 150 x 300 quadrat
trim = (200 - 150) / 2

i_bisections = c(13, 11, 9, 7, 5, 3)
n_quadrats = 2^i_bisections
domain = c(trim, 200 - trim, 0, 300) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$gx, dat$gy), n_quadrats, domain)

write.csv(comms, file='./data/ferp_comms.csv', row.names=FALSE)
