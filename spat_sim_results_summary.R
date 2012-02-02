## Purpose: to graphically summarize the simulated simical results
setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

files = dir('./sorensen')
shrtnames = c('bci','cocoli','sherman','serp')

fileNums = sapply(shrtnames,function(x)grep(x,files))
fileNums = lapply(fileNums,function(x) x[x %in% grep('C200',files)])

longnames = sub('sorensen_','',files[unlist(fileNums)])
longnames = sub('.Rdata','',longnames)
longnames = sub('_abu','',longnames)
longnames = sub('_binary','',longnames)
longnames = unique(longnames)

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

## now do one more round of averageing for cocoli and sherman
fixResults = function(results){
  out = vector('list',length=8)
  names(out) = names(results)[c(1:4,7:8,13:14)]
  out[1:2] = results[1:2]  ## bci 
  out[[3]] = rbind(results[[3]],results[[5]])  ## cocoli
  out[[4]] = rbind(results[[4]],results[[6]])  ## cocoli
  out[[5]] = rbind(results[[7]],results[[9]],results[[11]]) ## sherman
  out[[6]] = rbind(results[[8]],results[[10]],results[[12]]) ## sherman
  out[7:8] = results[13:14] ## serp
  return(out)
}

simSorBin = fixResults(simSorBinAvg)
simSorAbu = fixResults(simSorAbuAvg)
simVarBin = fixResults(simVarBinAvg)
simVarAbu = fixResults(simVarAbuAvg)

## last round of averaging
combine = vector('list',length=length(simSorBin))
combine[1:2] = NA
combine[[3]] = simSorBin[[3]]$Comm
combine[[4]] = simSorBin[[4]]$Comm
combine[[5]] = c(rep(1,142),rep(2,nrow(simSorBin[[5]])-142))
combine[[6]] = c(rep(1,142),rep(2,nrow(simSorBin[[6]])-142))
combine[7:8] = NA

simSorBinAvg = avgResults(simSorBin,combine)
simSorAbuAvg = avgResults(simSorAbu,combine)
simVarBinAvg = avgResults(simVarBin,combine)
simVarAbuAvg = avgResults(simVarAbu,combine)

## split results that were generated from log-series and those that were fixed abu
simSorBinLogSer = simSorBinAvg[c(1,3,5,7,9)]
simSorBinFixed = simSorBinAvg[c(2,4,6,8,10)]
simSorAbuLogSer = simSorAbuAvg[c(1,3,5,7,9)]
simSorAbuFixed = simSorAbuAvg[c(2,4,6,8,10)]
simVarBinLogSer = simVarBinAvg[c(1,3,5,7,9)]
simVarBinFixed = simVarBinAvg[c(2,4,6,8,10)]
simVarAbuLogSer = simVarAbuAvg[c(1,3,5,7,9)]
simVarAbuFixed = simVarAbuAvg[c(2,4,6,8,10)]

## Now bring in empirical Results
source('spat_empir_results_summary.R')

################ Abu Sorensen
par(mfrow=c(2,3))
results1 = empirSorAbuAvg
results2 = simSorAbuLogSer
results3 = simSorAbuFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,.4),main='bci')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,.05),main='cocoli')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,.06),main='sherman1')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(.001,.07),main='sherman2')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.1,.75),main='serp')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

################ Binary Sorensen
par(mfrow=c(2,3))
results1 = empirSorBinAvg
results2 = simSorBinLogSer
results3 = simSorBinFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,.4),main='bci')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,.06),main='cocoli')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,.07),main='sherman1')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(.001,.08),main='sherman2')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.1,1),main='serp')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

################ Abu varWithin
par(mfrow=c(2,3))
results1 = empirVarAbuAvg
results2 = simVarAbuLogSer
results3 = simVarAbuFixed
plot(Metric ~ Dist,data = results1[[1]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(35,2000),main='bci')
lines(Metric ~ Dist,data = results2[[1]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[1]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[2]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,1),main='cocoli')
lines(Metric ~ Dist,data = results2[[2]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[2]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(.001,4),main='sherman1')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[3]],subset= Comm == 13,col=1,lwd=2,type='l',ylim=c(.001,4),main='sherman2')
lines(Metric ~ Dist,data = results2[[3]],subset= Comm == 2,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[3]],subset= Comm == 2,col='dodgerblue',lwd=2,type='l')
#
plot(Metric ~ Dist,data = results1[[4]],subset= Comm == 1,col=1,lwd=2,type='l',ylim=c(1e3,6e4),main='serp')
lines(Metric ~ Dist,data = results2[[4]],subset= Comm == 1,col='red',lwd=2,type='l')
lines(Metric ~ Dist,data = results3[[4]],subset= Comm == 1,col='dodgerblue',lwd=2,type='l')

