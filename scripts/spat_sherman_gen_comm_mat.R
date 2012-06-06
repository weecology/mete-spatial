## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 3 of sherman
## which was surveyed in 1999 for each spatial scale of interest
## Metadata: sherman_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

## read in data from the 1999 census (i.e. census 3)

dat = read.csv('~/datasets/CTFSplots/sherman/sherman_census3_filtered.csv',
               header=TRUE)

uniSpeciesNames = as.character(sort(unique(dat$spcode)))
dat$spnum = match(dat$spcode,uniSpeciesNames)
S = max(dat$spnum) 

plot(dat$x,dat$y,pch='.')
abline(v=140,col='red')
abline(v=240,col='red')
abline(h=40,col='red')
abline(h=140,col='red')
abline(h=440,col='red')

## split plot into three quadrats (two 200 x 100 m rectangles and one 140 x 140m)
dat1 = dat[dat$x >= 140 & dat$y >= 240, ]  ## one of the 200 x 100 m rectangles
dat2 = dat[dat$x >= 140 & dat$y < 240, ] ## the other 200 x 100 m rectangle
dat3 = dat[dat$x < 140, ]  ## the 140 x 140 square

plot(dat1$x,dat1$y,col='red',pch='.',ylim=range(dat$y),xlim=range(dat$x))
points(dat2$x,dat2$y,col='blue',pch='.')
points(dat3$x,dat3$y,col='green3',pch='.')

i_bisections = c(13, 11, 9, 7, 5, 3)
n_quadrats = 2^i_bisections
domain1 = c(140,240,240,440) # spatial domain in meters defined here
domain2 = c(140,240,40,240)

## generate a site x species matrix for each spatial scale
comms1 = make_comm_matrix(dat1$spnum, S, cbind(dat1$x, dat1$y), n_quadrats,
                          domain1)
comms2 = make_comm_matrix(dat2$spnum, S, cbind(dat2$x, dat2$y), n_quadrats,
                          domain2)

## now for dat3 the square plot
i_bisections = c(12, 10, 8, 6, 4)
n_quadrats = 2^i_bisections
domain3 = c(0,140,0,140)

comms3 = make_comm_matrix(dat3$spnum, S, cbind(dat3$x, dat3$y), n_quadrats,
                          domain3)

## output results
write.csv(comms1,file='./data/sherman1_comms.csv',row.names=F)
write.csv(comms2,file='./data/sherman2_comms.csv',row.names=F)
write.csv(comms3,file='./data/sherman3_comms.csv',row.names=F)
