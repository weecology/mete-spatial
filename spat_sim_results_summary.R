## Purpose: to graphically summarize the simulated simical results
setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

files = dir('./sorensen')
shrtnames = c('bci','cocoli1','cocoli2','cross','serp','sherman1','sherman2','sherman3')

fileNums = sapply(shrtnames,function(x)grep(x,files))
fileNums = lapply(fileNums,function(x) x[x %in% grep('C200',files)])

longnames = sub('sorensen_','',files[unlist(fileNums)])
longnames = sub('.Rdata','',longnames)
longnames = sub('_abu','',longnames)
longnames = sub('_binary','',longnames)
longnames = unique(longnames)

## fix problem with no having a crosstimbers logseries fit
longnames = c(longnames[1:7],longnames[7:15])

simBin = getResults(longnames,'sorensen','binary')
simAbu = getResults(longnames,'sorensen','abu')
simSorBin = reshapeResults(simBin,'sorensen')
simSorAbu = reshapeResults(simAbu,'sorensen')

simBin = getResults(longnames,'varWithin','binary')
simAbu = getResults(longnames,'varWithin','abu')
simVarBin = reshapeResults(simBin,'varWithin')
simVarAbu = reshapeResults(simAbu,'varWithin')

combine = vector('list',length=length(simSorBin))
for(i in seq_along(combine))
  combine[[i]] = rep(1,nrow(simSorBin[[i]]))

simSorBinAvg = avgResults(simSorBin,combine)
simSorAbuAvg = avgResults(simSorAbu,combine)
simVarBinAvg = avgResults(simVarBin,combine)
simVarAbuAvg = avgResults(simVarAbu,combine)

## split results that were generated from log-series and those that were fixed abu
simSorBinLogSer = simSorBinAvg[c(1,3,5,7,9,11,13,15)]
simSorBinFixed = simSorBinAvg[c(2,4,6,8,10,12,14,16)]
simSorAbuLogSer = simSorAbuAvg[c(1,3,5,7,9,11,13,15)]
simSorAbuFixed = simSorAbuAvg[c(2,4,6,8,10,12,14,16)]
simVarBinLogSer = simVarBinAvg[c(1,3,5,7,9,11,13,15)]
simVarBinFixed = simVarBinAvg[c(2,4,6,8,10,12,14,16)]
simVarAbuLogSer = simVarAbuAvg[c(1,3,5,7,9,11,13,15)]
simVarAbuFixed = simVarAbuAvg[c(2,4,6,8,10,12,14,16)]

## Now bring in empirical Results
source('spat_empir_results_summary.R')

################ Abu Sorensen
pdf('empir_sim_sor_arith.pdf')
par(mfrow=c(2,4))
results1 = empirSorAbu
results2 = simSorAbuLogSer
results3 = simSorAbuFixed
ylims = rbind(c(0.001,.4),c(.001,.05),c(.001,.05),c(.001,.4),
              c(.001,.7),c(.001,.05),c(.001,.05),c(.001,.05))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l',
       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}

################ Binary Sorensen
par(mfrow=c(2,4))
results1 = empirSorBin
results2 = simSorBinLogSer
results3 = simSorBinFixed
ylims = rbind(c(0.001,.4),c(.001,.05),c(.001,.05),c(.001,.6),
              c(.3,1),c(.001,.05),c(.001,.05),c(.001,.05))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l',
       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}
dev.off()
################ Abu varWithin
par(mfrow=c(2,4))
results1 = empirVarAbu
results2 = simVarAbuLogSer
results3 = simVarAbuFixed
ylims = rbind(c(35,2000),c(.001,.05),c(.001,.05),c(.001,.4),
              c(.3,.7),c(.001,.05),c(.001,.05),c(.001,.05))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l')
#       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}

################ binary varWithin
par(mfrow=c(2,4))
results1 = empirVarBin
results2 = simVarBinLogSer
results3 = simVarBinFixed
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='l')
#       ylim=ylims[i,],main=shrtnames[i])
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='dodgerblue',
        lwd=2,type='l')
}

dev.off()
############################################ log -log
################ Abu Sorensen
pdf('empir_sim_sor_loglog.pdf')
par(mfrow=c(2,4))
results1 = empirSorAbu
results2 = simSorAbuLogSer
results3 = simSorAbuFixed
ylims = rbind(c(0.01,1),c(.001,1),c(.001,1),c(.001,.4),
              c(.01,1),c(.003,1),c(.003,1),c(.003,1))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='n',
       ylim=ylims[i,],main=shrtnames[i],log='xy')
  tmp = results2[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='palegreen')
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='green3',lwd=2,
        type='l')
  tmp = results3[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='dodgerblue')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='blue',
        lwd=2,type='l')
  lines(Exp ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  points(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,pch=19)
}

################ Binary Sorensen
par(mfrow=c(2,4))
results1 = empirSorBin
results2 = simSorBinLogSer
results3 = simSorBinFixed
ylims = rbind(c(0.01,1),c(.001,1),c(.001,1),c(.001,.4),
              c(.01,1),c(.003,1),c(.003,1),c(.003,1))
for(i in seq_along(shrtnames)){
  plot(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,lwd=2,type='n',
       ylim=ylims[i,],main=shrtnames[i],log='xy')
  tmp = results2[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='palegreen')
  lines(Metric ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='green3',lwd=2,
        type='l')
  tmp = results3[[i]]
  polygon(c(tmp$Dist,rev(tmp$Dist)),c(tmp$MetricHi,rev(tmp$MetricLo)),border=NA,
          col='dodgerblue')
  lines(Metric ~ Dist,data = results3[[i]],subset= Comm == 1|Comm == 10,col='blue',
        lwd=2,type='l')
  lines(Exp ~ Dist,data = results2[[i]],subset= Comm == 1|Comm == 10,col='red',lwd=2,
        type='l')
  points(Metric ~ Dist,data = results1[[i]],subset= Comm == 1 |Comm == 10,col=1,pch=19)
}
dev.off()

par(mfrow=c(1,1))
plot(1:10,1:10,frame.plot=F,axes=F,xlab='',ylab='',type='n')
legend('center',c('Empirical','METE','METE fixed abu','RP'),col=c('black','green3',
       'blue','red'),cex=2,pch=c(19,rep(NA,3)),lty=c(NA,rep(1,3)),lwd=c(NA,rep(2,4)),bty='n')



################ Abu varWithin
par(mfrow=c(2,3))
results1 = empirVarAbuAvg
results2 = simVarAbuLogSer
results3 = simVarAbuFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(35,2000),main='bci',log='xy')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.1,1),main='cocoli',log='xy')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.1,4),main='sherman1',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(1,4),main='sherman2',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(1e3,6e4),main='serp',log='xy')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

################ binary varWithin
par(mfrow=c(2,3))
results1 = empirVarBinAvg
results2 = simVarBinLogSer
results3 = simVarBinFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(4,13),main='bci',log='xy')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.35,.55),main='cocoli',log='xy')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.6,1),main='sherman1',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(1,1.6),main='sherman2',log='xy')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(1,4),main='serp',log='xy')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

dev.off()
