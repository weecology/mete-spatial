
setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

dat = read.csv('./data/filtered_data/luqu_census4_filtered.csv')

uniSpeciesNames = as.character(sort(unique(dat$SPCODE)))
dat$spnum = match(dat$SPCODE, uniSpeciesNames)
S = max(dat$spnum) 

range(dat$GX) ## 320 max
range(dat$GY) ## 500 max
## change into a single 250 x 500 quadrat
trim = (320 - 250) / 2

i_bisections = c(13, 11, 9, 7, 5, 3)
n_quadrats = 2^i_bisections
domain = c(trim, 320 - trim, 0, 500) # spatial domain in meters defined here

## generate a site x species matrix for each spatial scale

comms = make_comm_matrix(dat$spnum, S, cbind(dat$GX, dat$GY), n_quadrats, domain)

write.csv(comms, file='./data/luquillo_comms.csv', row.names=FALSE)
