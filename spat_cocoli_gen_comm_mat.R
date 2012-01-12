## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 3 of cocoli
## which was surveyed in 1998 for each spatial scale of interest
## Metadata: cocoli_files.doc 
## There is a visual diagram of the Sherman plot in the following publication,
## which is fairly close to the cocoli plot except for two 40 m strips 
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Niño dry season. Journal of Tropical Ecology, 20: 51-72.

## read in data from the 1998 census (i.e. census 3)

dat = read.csv('/home/danmcglinn/CTFSplots/cocoli/cocoli_census3_filtered.csv',header=TRUE)

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

len = function(Nbisect){
  2^floor(Nbisect/2) 
}

wid = function(Nbisect){
  sapply(Nbisect,function(Nbisect){
  if(Nbisect %% 2 == 0)
   len(Nbisect)/2
  else
   len(Nbisect)
  })
}

## for square quadrats the lengths are in meters defined here
quadLen = round(shortSide/ wid(c(14,12,10,8,6,4)), 2)
quadN = round(shortSide/quadLen) * round(longSide/quadLen)
## generate a site x species matrix for each spatial scale

comms = matrix(NA, nrow=sum(quadN)*2 ,ncol=S+3)
colnames(comms) = c('grain','x','y',paste('sp',1:S,sep=''))
irow = 1

## work with dat1 first 
minx = 0
maxx = 100
miny = 100
maxy = 300

for(A in seq_along(quadLen)){
  xbreaks = seq(minx,maxx,quadLen[A])
  ybreaks = seq(miny,maxy,quadLen[A]) 
  for (x in 1:(length(xbreaks)-1)) {
    for (y in 1:(length(ybreaks)-1)) {
      inQuad =  xbreaks[x] <= dat1$x & dat1$x < xbreaks[x+1] & 
                ybreaks[y] <= dat1$y & dat1$y < ybreaks[y+1]   
      comms[irow, c(1:3)] = c(paste(round(quadLen[A]^2),'_1',sep=''),x,y)
      comms[irow, -c(1:3)] = as.integer(table(c(dat1$spnum[inQuad],1:S)) - 1)
      irow = irow + 1 
    }
  }
}

## work with dat2 now
minx = 0
maxx = 200
miny = 0
maxy = 100

for(A in seq_along(quadLen)){
  xbreaks = seq(minx,maxx,quadLen[A])
  ybreaks = seq(miny,maxy,quadLen[A]) 
  for (x in 1:(length(xbreaks)-1)) {
    for (y in 1:(length(ybreaks)-1)) {
      inQuad =  xbreaks[x] <= dat2$x & dat2$x < xbreaks[x+1] & 
                ybreaks[y] <= dat2$y & dat2$y < ybreaks[y+1]   
      comms[irow, c(1:3)] = c(paste(round(quadLen[A]^2),'_2',sep=''),x,y)
      comms[irow, -c(1:3)] = as.integer(table(c(dat2$spnum[inQuad],1:S)) - 1)
      irow = irow + 1 
    }
  }
}

write.csv(comms,file='/home/danmcglinn/maxent/spat/data/cocoli_comms.csv',
          row.names=FALSE)
save(comms,file='/home/danmcglinn/CTFSplots/cocoli/cocoli_comms_census3.Rdata')

     