## Purpose: to graphically summarize the empirical results
setwd('/home/danmcglinn/maxent/spat')
source('spat_sim_vario_func.R')

shrtnames = c('bci','cocoli1','cocoli2','cross','serp','sherman1','sherman2','sherman3')

empirBin = getResults(shrtnames,'sorensen','binary')
empirAbu = getResults(shrtnames,'sorensen','abu')
empirSorBin = reshapeResults(empirBin,'sorensen')
empirSorAbu = reshapeResults(empirAbu,'sorensen')

empirBin = getResults(shrtnames,'varWithin','binary')
empirAbu = getResults(shrtnames,'varWithin','abu')
empirVarBin = reshapeResults(empirBin,'varWithin')
empirVarAbu = reshapeResults(empirAbu,'varWithin')

## Average cocoli1 & 2 and sherman 1 & 2
#cocoliSorBin = (empirSorBin[[2]] + empirSorBin[[3]])/2
#shermanSorBin = (empirSorBin[[6]] + empirSorBin[[7]])/2
#empirSorBin = empirSorBin[-c(3,7)]
#empirSorBin[[2]] = cocoliSorBin
#empirSorBin[[6]] = shermanSorBin
#names(empirSorBin) = c('bci','cocoli','cross','serp','sherman','sherman3')


#pdf('spat_empir_sor_bin_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirSorBin)
plotEmpir(empirSorBin,log='xy')
#dev.off()
#pdf('spat_empir_sor_abu_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirSorAbu)
plotEmpir(empirSorAbu,log='xy')
#dev.off()
#pdf('spat_empir_varWithin_bin_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirVarBin)
plotEmpir(empirVarBin,log='xy')
#dev.off()
#pdf('spat_empir_varWithin_abu_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirVarAbu)
plotEmpir(empirVarAbu,log='xy')
#dev.off()

#vGridExpAvg = apply(vGridExp,1,mean)
#vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))


