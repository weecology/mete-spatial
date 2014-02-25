## there are two methods to how a complete spatial random (CSR)
## expectation can be generated
## one can randomize samples or one can randomize individuals
## they will produce different expected amounts of species turnover
## if observed individuals are not distributed according to a homogenous
## Poisson process.

library(vegan)

source('./scripts/spat_functions.R')

site = 'ferp'
dat = read.csv(paste('./data/', site, '_comms.csv', sep=''))

unique(dat[,1])

dat_sub = dat[dat$grain == 87.89, ]
comm = as.matrix(dat_sub[ , -(1:3)])
coords = dat_sub[ , 2:3]

dim(dat_sub)
log2(nrow(comm))

csiml = 1 - vegdist(comm)
gdist = dist(coords)

v = vario_bisect(comm, coords, distance.metric='bray')

nperm = 20
vperm = matrix(NA, ncol=nperm, nrow=length(v$vario$exp.var))

for(i in 1:nperm) {
  ## generate new SSAD
  nquad = nrow(comm)
  comm_tmp = comm
  for (j in 1:ncol(comm)) {
    rand_samp = sample(nquad, size=sum(comm[ , j]), replace=T)
    comm_tmp[ , j] = as.numeric(table(c(rand_samp, 1:nquad)) - 1)
  }
  vperm[ , i] = vario_bisect(comm_tmp, coords, distance.metric='bray')$vario$exp.var
}

vperm_avg = apply(vperm, 1, mean)

par(mfrow=c(1,2))
plot(v$vario$Dist, 1 - v$vario$exp.var, type='o', log='xy', ylim=c(.2, 0.7))
lines(v$vario$Dist, 1 - vperm2_avg, lwd=3, col='red')
abline(h=mean(csiml, na.rm=T), col='blue', lwd=3)
plot(1:10, 1:10, type='n', axes=F, xlab='', ylab='')
legend('center', c('observed', 'Sample-based RPM', 'Individual-based RPM'),
       col=c('black','blue','red'), lty=1, pch=c(1, NA, NA), cex=.55, bty='n')



