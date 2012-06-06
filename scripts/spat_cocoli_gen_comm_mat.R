## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 3 of cocoli
## which was surveyed in 1998 for each spatial scale of interest
## Metadata: cocoli_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## which is fairly close to the Cocoli plot except for two 40 m strips 
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

## read in data from the 1998 census (i.e. census 3)

dat = read.csv('~/datasets/CTFSplots/cocoli/cocoli_census3_filtered.csv',
               header=TRUE)

uniSpeciesNames = as.character(sort(unique(dat$spcode)))
dat$spnum = match(dat$spcode,uniSpeciesNames)
S = max(dat$spnum) 

plot(dat$x,dat$y,pch='.')

## split plot into two 200 x 100 m rectangles
dat1 = dat[dat$y>=100,]
dat2 = dat[dat$y<100,]

plot(dat1$x,dat1$y,col='red',pch='.',ylim=range(dat$y),xlim=range(dat$x))
points(dat2$x,dat2$y,col='blue',pch='.')

i_bisections = c(13, 11, 9, 7, 5, 3)
n_quadrats = 2^i_bisections
domain1 = c(0,100,100,300) # spatial domain in meters defined here
domain2 = c(0,200,0,100) 

## generate a site x species matrix for each spatial scale

comms1 = make_comm_matrix(dat$spnum, S, cbind(dat$x, dat$y), n_quadrats,
                          domain1)
comms2 = make_comm_matrix(dat$spnum, S, cbind(dat$x, dat$y), n_quadrats,
                          domain2)

## output results
write.csv(comms1,file='./data/cocoli1_comms.csv',row.names=F)
write.csv(comms2,file='./data/cocoli2_comms.csv',row.names=F)
