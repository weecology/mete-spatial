## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 3 of sherman
## which was surveyed in 1999 for each spatial scale of interest
## Metadata: sherman_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

source('/home/danmcglinn/maxent/spat/spat_sim_vario_func.R')

## read in data from the 1999 census (i.e. census 3)

dat = read.csv('/home/danmcglinn/CTFSplots/sherman/sherman_census3_filtered.csv',
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
dat1 = dat[dat$x>=140 & dat$y>=240,]  ## one of the 200 x 100 m rectangles
dat2 = dat[dat$x>=140 & dat$y<240,] ## the other 200 x 100 m rectangle
dat3 = dat[dat$x<140,]  ## the 140 x 140 square

plot(dat1$x,dat1$y,col='red',pch='.',ylim=range(dat$y),xlim=range(dat$x))
points(dat2$x,dat2$y,col='blue',pch='.')
points(dat3$x,dat3$y,col='green3',pch='.')

## work with dat1 first 

shortSide = 100
longSide = 200
## for square quadrats the lengths are in meters defined here
nPixels = wid(c(14,12,10,8,6,4))
quadLen = shortSide/ nPixels
quadN = nPixels^2 * (longSide / shortSide)
## generate a site x species matrix for each spatial scale
comms1 = makeCommMat(dat1$spnum,S,cbind(dat1$x,dat1$y),quadLen,quadN,
                     c(140,240,240,440),'_1')

## work with dat2 now

shortSide = 100
longSide = 200
## for square quadrats the lengths are in meters defined here
nPixels = wid(c(14,12,10,8,6,4))
quadLen = shortSide/ nPixels
quadN = nPixels^2 * (longSide / shortSide)
## generate a site x species matrix for each spatial scale
comms2 = makeCommMat(dat2$spnum,S,cbind(dat2$x,dat2$y),quadLen,quadN,
                     c(140,240,40,240),'_2')

## work with dat3 now

shortSide = 140
longSide = 140
## for square quadrats the lengths are in meters defined here
nPixels = wid(c(14,12,10,8,6,4))
quadLen = shortSide/ nPixels
quadN = nPixels^2 * (longSide / shortSide)
## generate a site x species matrix for each spatial scale
comms3 = makeCommMat(dat3$spnum,S,cbind(dat3$x,dat3$y),quadLen,quadN,
                     c(0,140,0,140),'_3')

comms = rbind(comms1,comms2,comms3)

write.csv(comms,file='/home/danmcglinn/maxent/spat/data/sherman_comms.csv',
          row.names=FALSE)
save(comms,file='/home/danmcglinn/CTFSplots/sherman/sherman_comms_census3.Rdata')

