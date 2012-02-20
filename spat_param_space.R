setwd('c:/users/dan mcglinn/documents/lab data/maxent/spat/')
source('spat_sim_vario_func.R')

files = dir('./sorensen/')

S <- round(10^seq(log10(10),log10(100),length.out=20))
N <- round(10^seq(log10(120),log10(5e5),length.out=20))

results = vector('list',length(S)*length(N))
icount = 0
for(s in S){
  for(n in N){
    icount = icount + 1
    fileSuffix = paste('_S',s,'_N',n,sep='')
    fileName = files[grep(fileSuffix,files)]
    if(length(fileName) == 0){
      next 
    }
    load(file.path('./sorensen/',fileName))
    results[[icount]] = sorensen
    rm(sorensen)
  }
}

simSorAbu = reshapeResults(results,'sorensen')

combine = vector('list',length=length(simSorAbu))
for(i in seq_along(combine)){
  if(!is.null(simSorAbu[[i]]))
    combine[[i]] = rep(1,nrow(simSorAbu[[i]]))
}

simSorAbuAvg = avgResults(simSorAbu,combine)

results = simSorAbuAvg
stats = array(NA,dim=c(3,length(S),length(N)))
dimnames(stats)[[1]] = c('b0','b1','r2')
dimnames(stats)[[2]] = S
dimnames(stats)[[3]] = N
i = 0
for(s in seq_along(S)){
  for(n in seq_along(N)){
    i = i + 1
    if(is.null(results[[i]]))
      next
    mod = lm(log10(Metric) ~ log10(Dist),data = results[[i]])
    stats[1:2,s,n] = coef(mod)
    stats[3,s,n] = summary(mod)$r.squ
  }
}

par(mfrow=c(1,3))
image(log10(S),log10(N[1:11]),stats[1,,1:11],xlab='',ylab='',bg='red')
image(log10(S),log10(N[1:11]),stats[2,,1:11],xlab='',ylab='')
image(log10(S),log10(N[1:11]),stats[3,,1:11],xlab='',ylab='')

par(mfrow=c(1,3))
for(i in 1:3){
  tmp = as.vector(stats[i,,1:11])
  x = round(seq(min(tmp,na.rm=T),max(tmp,na.rm=T),length.out=7),3)
  image(matrix(x),axes=F)
  axis(side=1,at=seq(0,1,length.out=7),labels=x,tick=F,cex.axis=2)
}

svals = rep(S,each = length(N))
nvals = rep(N,length(S))
ratio = log10(nvals/svals)
ratio = array(ratio,dim=c(length(S),length(N)))

pdf('dd_coef_n_over_s.pdf')
lwd = 3
par(mfrow=c(2,3))
xlab = 'log10(N/S)'
plot(ratio,stats[1,,],xlab=xlab,ylab='Intercept',pch=19)
plot(ratio,stats[1,,],xlab=xlab,ylab='Intercept',type='n')
for(s in seq_along(S))
  lines(ratio[s,],stats[1,s,],col='palevioletred',lwd=lwd)
plot(ratio,stats[1,,],xlab=xlab,ylab='Intercept',type='n')
for(n in seq_along(N))
  lines(ratio[,n],stats[1,,n],col='dodgerblue',lwd=lwd)
##
plot(ratio,stats[2,,],xlab=xlab,ylab='Slope',pch=19)
plot(ratio,stats[2,,],xlab=xlab,ylab='Slope',type='n')
for(s in seq_along(S))
  lines(ratio[s,],stats[2,s,],col='palevioletred',lwd=lwd)
plot(ratio,stats[2,,],xlab=xlab,ylab='Slope',type='n')
for(n in seq_along(N))
  lines(ratio[,n],stats[2,,n],col='dodgerblue',lwd=lwd)
dev.off()
#### same as above but for presentation
lwd = 3
par(mfrow=c(2,3))
plot(ratio,stats[1,,],xlab='',ylab='',pch=19)
plot(ratio,stats[1,,],xlab='',ylab='',type='n')
for(s in seq_along(S))
  lines(ratio[s,],stats[1,s,],col='palevioletred',lwd=lwd)
plot(ratio,stats[1,,],xlab='',ylab='',type='n')
for(n in seq_along(N))
  lines(ratio[,n],stats[1,,n],col='dodgerblue',lwd=lwd)
##
plot(ratio,stats[2,,],xlab='',ylab='',pch=19)
plot(ratio,stats[2,,],xlab='',ylab='',type='n')
for(s in seq_along(S))
  lines(ratio[s,],stats[2,s,],col='palevioletred',lwd=lwd)
plot(ratio,stats[2,,],xlab='',ylab='',type='n')
for(n in seq_along(N))
  lines(ratio[,n],stats[2,,n],col='dodgerblue',lwd=lwd)


###Scale collapse attempt
ratio = log10(nvals/svals) - (svals *log10(svals))
ratio = array(ratio,dim=c(length(S),length(N)))
xlab = 'log10(N/S) - (S * log10(S))'
par(mfrow=c(1,3))
plot(ratio,stats[2,,],xlab=xlab,ylab='Slope',pch=19,xlim=c(-60,0))
plot(ratio,stats[2,,],xlab=xlab,ylab='Slope',type='n',xlim=c(-60,0))
for(s in seq_along(S))
  lines(ratio[s,],stats[2,s,],col='palevioletred',lwd=lwd)
plot(ratio,stats[2,,],xlab=xlab,ylab='Slope',type='n',xlim=c(-60,0))
for(n in seq_along(N))
  lines(ratio[,n],stats[2,,n],col='dodgerblue',lwd=lwd)





