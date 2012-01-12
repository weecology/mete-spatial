## Author: Dan McGlinn
## Purpose: to create site x species matrices for census 7
## which was surveyed in 2010 for each spatial scale of interest
## Metadata: bci50ha.doc 

## read in data from the 2010 census (i.e. census 7)

load('/home/danmcglinn/CTFSplots/BCI/bci_census7.Rdata')

uniSpeciesNames = as.character(sort(unique(dat$Latin)))

dat$spnum = match(dat$Latin,uniSpeciesNames)

S = max(dat$spnum) 

attach(dat)

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
quadLen = round(500 / wid(c(14,12,10,8,6,4)), 2)

quadN = round(1e3/quadLen) * round(5e2/quadLen)

## generate a site x species matrix for each spatial scale

comms = matrix(NA, nrow=sum(quadN) ,ncol=S+3)
colnames(comms) = c('grain','x','y',paste('sp',1:S,sep=''))
irow = 1
     
for(A in seq_along(quadLen)){
  xbreaks <- seq(0,1000,quadLen[A])
  ybreaks <- seq(0,500,quadLen[A]) 
  for (x in 1:(length(xbreaks)-1)) {
    for (y in 1:(length(ybreaks)-1)) {
      inQuad =  xbreaks[x] <= gx & gx < xbreaks[x+1] & 
                ybreaks[y] <= gy & gy < ybreaks[y+1]   
      comms[irow, c(1:3)] = c(round(quadLen[A]^2),x,y)
      comms[irow, -c(1:3)] = as.integer(table(c(spnum[inQuad],1:S)) - 1)
      irow = irow + 1 
    }
  }
}

write.csv(comms,file='/home/danmcglinn/maxent/spat/data/bci_comms.csv',
          row.names=FALSE)
save(comms,file='/home/danmcglinn/CTFSplots/BCI/bci_comms_census7.Rdata')

     