## Purpose: to graphically summarize the empirical results
setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

shrtnames = c('bci','cocoli','sherman','cross','serp')

empirBin = getResults(shrtnames,'sorensen','binary')
empirAbu = getResults(shrtnames,'sorensen','abu')
empirSorBin = reshapeResults(empirBin,'sorensen')
empirSorAbu = reshapeResults(empirAbu,'sorensen')

empirBin = getResults(shrtnames,'varWithin','binary')
empirAbu = getResults(shrtnames,'varWithin','abu')
empirVarBin = reshapeResults(empirBin,'varWithin')
empirVarAbu = reshapeResults(empirAbu,'varWithin')

combine = vector('list',length=length(empirSorBin))
combine[[1]] = NA
combine[[2]] = ifelse(empirSorBin$cocoli$Comm > 6,empirSorBin$cocoli$Comm - 6,
                      empirSorBin$cocoli$Comm)
combine[[3]] = ifelse(empirSorBin$sherman$Comm > 6 & empirSorBin$sherman$Comm < 13,
                      empirSorBin$sherman$Comm - 6, empirSorBin$sherman$Comm)
combine[[4]] = NA

empirSorBinAvg = avgResults(empirSorBin,combine)
empirSorAbuAvg = avgResults(empirSorAbu,combine)
empirVarBinAvg = avgResults(empirVarBin,combine)
empirVarAbuAvg = avgResults(empirVarAbu,combine)

#pdf('spat_empir_sor_bin_curves.pdf')
par(mfrow=c(2,2))
plotEmpir(empirSorBinAvg)
plotEmpir(empirSorBinAvg,log='xy')
#dev.off()
#pdf('spat_empir_sor_abu_curves.pdf')
par(mfrow=c(2,2))
plotEmpir(empirSorAbuAvg)
plotEmpir(empirSorAbuAvg,log='xy')
#dev.off()
#pdf('spat_empir_varWithin_bin_curves.pdf')
par(mfrow=c(2,2))
plotEmpir(empirVarBinAvg)
plotEmpir(empirVarBinAvg,log='xy')
#dev.off()
#pdf('spat_empir_varWithin_abu_curves.pdf')
par(mfrow=c(2,2))
plotEmpir(empirVarAbuAvg)
plotEmpir(empirVarAbuAvg,log='xy')
#dev.off()

#vGridExpAvg = apply(vGridExp,1,mean)
#vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))


