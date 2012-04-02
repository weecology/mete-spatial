
setwd('./maxent/spat/')
library(spatstat)

dat = read.csv('C:/Users/Dan McGlinn/Dropbox/Public/cross1998_filtered.csv')
## drop points that have the same x and y coordinates
dim(dat)
dat.uni = unique(data.frame(sp = dat$sp, gx = dat$x, gy = dat$y))
dim(dat.uni)
## drop low abundance species
dat.uni = dat.uni[dat.uni$sp == 'QUMA'| dat.uni$sp=='QUST'| dat.uni$sp=='CECO',]
dat.uni$sp = factor(dat.uni$sp,levels=c('QUMA','QUST'))

head(dat)
?ppp

dat.ppp = ppp(dat.uni$gx,dat.uni$gy,window=owin(c(0,200),c(0,200),unitname=c('meter','meters')),
              marks=data.frame(sp=dat.uni$sp))

summary(dat.ppp)

plot(Kest(dat.ppp))

qust.ppp = subset(dat.ppp,dat.ppp$marks$sp == 'QUST')
quma.ppp = subset(dat.ppp,dat.ppp$marks$sp == 'QUMA')
K = Kest(qust.ppp)
par(mfrow=c(1,2))
plot(Kest(qust.ppp))
plot(Kest(quma.ppp))

K.qust = Kmulti(dat.ppp,dat.ppp$marks=='QUST',dat.ppp$marks=='QUST')
K.quma = Kmulti(dat.ppp,dat.ppp$marks=='QUMA',dat.ppp$marks=='QUMA')

par(mfrow=c(1,2))
plot(K.qust,legend=F)
plot(K.quma,legend=F)

source('Shen et al. (2009) - rcode.R')
n = 100
dat.tst = data.frame(sp = rep('a',n),gx = runif(n)*200, gy = runif(n)*200)
dat.tst = rbind(dat.tst,
                data.frame(sp = rep('b',n),gx = runif(n)*200, gy = runif(n)*200))

SAR.Obs = SAR(dat.tst,10,c(200,200),process='Real')
SAR.Poi = SAR(dat.tst,10,c(200,200),process='Hom.Po')
SAR.Th = SAR(dat.tst,10,c(200,200),process='Hom.Th')

plot(SAR.Obs$A,SAR.Obs$S.mean,log='xy',type='o',pch=19)
points(SAR.Poi$A,SAR.Poi$S.mean,pch=1,col='blue',lwd=2)
points(SAR.Th$A,SAR.Th$S.mean,pch=1,col='red',lwd=2)

### reformat crosstimbers to work for this 
cross = data.frame(sp = dat.uni$sp,gx = dat.uni$gx, gy = dat.uni$gy)

SAR.Obs = SAR(cross,10,c(200,200),process='Real')
SAR.Poi = SAR(cross,10,c(200,200),process='Hom.Po')
SAR.Th = SAR(cross,10,c(200,200),process='Hom.Th')




