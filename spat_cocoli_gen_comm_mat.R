## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 3 of cocoli
## which was surveyed in 1998 for each spatial scale of interest
## Metadata: cocoli_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## which is fairly close to the Cocoli plot except for two 40 m strips 
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

source('/home/danmcglinn/maxent/spat/spat_sim_vario_func.R')

## read in data from the 1998 census (i.e. census 3)

dat = read.csv('/home/danmcglinn/CTFSplots/cocoli/cocoli_census3_filtered.csv',
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

shortSide = 100
longSide = 200

## for square quadrats the lengths are in meters defined here
nPixels = wid(c(14,12,10,8,6,4))
quadLen = shortSide/ nPixels
quadN = nPixels^2 * (longSide / shortSide)
## generate a site x species matrix for each spatial scale

## work with dat1 first 
comms1 = makeCommMat(dat1$spnum,S,cbind(dat1$x,dat1$y),quadLen,quadN,
                     c(0,100,100,300),'_1')

## work with dat2 now
comms2 = makeCommMat(dat2$spnum,S,cbind(dat2$x,dat2$y),quadLen,quadN,
                     c(0,200,0,100),'_2')

comms = rbind(comms1,comms2)

write.csv(comms,file='/home/danmcglinn/maxent/spat/data/cocoli_comms.csv',
          row.names=FALSE)
save(comms,file='/home/danmcglinn/CTFSplots/cocoli/cocoli_comms_census3.Rdata')

