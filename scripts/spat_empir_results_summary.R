## Purpose: to graphically summarize the empirical results
setwd('./maxent/spat')
source('./scripts/spat_sim_vario_func.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

merge_drop = function(results){
  result_names = names(results)
  to_drop = match(c('cocoli2','sherman2','sherman3'), result_names)
  cocoli_results = (results$cocoli1 + results$cocoli2)/2
  sherman_results = (results$sherman1 + results$sherman2)/2
  results = results[-to_drop]
  results[['cocoli1']] = cocoli_results
  results[['sherman1']] = sherman_results
  return(results)
}

empirBin = getResults(shrtnames,'sorensen','binary')
empirAbu = getResults(shrtnames,'sorensen','abu')
empirSorBin = reshapeResults(empirBin,'sorensen')
empirSorAbu = reshapeResults(empirAbu,'sorensen')
## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirSorBin = merge_drop(empirSorBin)
empirSorAbu = merge_drop(empirSorAbu)

par(mfrow=c(2,4))
plotEmpir(empirSorAbu,log='x', type='o')
plotEmpir(empirSorAbu,log='x', type='o',quants=T)

empirBin = getResults(shrtnames,'varWithin','binary')
empirAbu = getResults(shrtnames,'varWithin','abu')
empirVarBin = reshapeResults(empirBin,'varWithin')
empirVarAbu = reshapeResults(empirAbu,'varWithin')
## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirVarBin = merge_drop(empirVarBin)
empirVarAbu = merge_drop(empirVarAbu)

par(mfrow=c(2,4))
plotEmpir(empirVarAbu,log='x', type='o')
plotEmpir(empirVarAbu,log='x', type='o',quants=T)


##
commName = 'bryan'
par(mfrow=c(1,2))
plotEmpir(empirSorAbu[commName],log='x', type='o', quants=T)
plotEmpir(empirVarAbu[commName],log='x', type='o', quants=T)


#pdf('spat_empir_sor_bin_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirSorBin, quants=TRUE, type='o')
plotEmpir(empirSorBin, quants=FALSE, type='o')

plotEmpir(empirSorBin[3], log='x', quants=T)
#dev.off()
#pdf('spat_empir_sor_abu_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirSorAbu, quants=TRUE)
plotEmpir(empirSorAbu, quants=FALSE)
plotEmpir(empirSorAbu,log='x', type='o')
#dev.off()
#pdf('spat_empir_varWithin_bin_curves.pdf')
par(mfrow=c(3,4))
plotEmpir(empirVarBin, quants=TRUE)
plotEmpir(empirVarBin,log='xy')
#dev.off()
#pdf('spat_empir_varWithin_abu_curves.pdf')
par(mfrow=c(2,4))
plotEmpir(empirVarAbu)
plotEmpir(empirVarAbu,log='xy')
#dev.off()

#vGridExpAvg = apply(vGridExp,1,mean)
#vGridExpQt = apply(vGridExp,1,function(x)quantile(x,c(.025,.975)))

## finding the optimal bin sizes
## bci
comm = 'bci'
results = empirSorAbu[comm]
unigrains = unique(results[[1]]$Comm)
col = palette()[-1]
par(mfrow=c(1,1))
plotEmpir(results, type='o', log='xy')
for (i in seq_along(unigrains))
  abline(v = 2 * sqrt(unigrains[i]),col=col[i])
empirSorAbu[comm]
####
r2 = results[[1]][results[[1]]$Comm==unigrains[1],]
breaks = exp(seq(0, log(max(r2$Dist) * sqrt(unigrains[1])),
                        length.out=20))
H = r2$Dist
for (i in 1:(length(breaks) - 1))
  H[H >= breaks[i] & H < breaks[i + 1]] = breaks[i]

en.split = split(r2$Metric.50 * r2$N , H)
n.split = split(r2$N, H)
dn.split = split(r2$Dist * sqrt(unigrains[1]) * r2$N, H)
e.mean = sapply(en.split, sum) / sapply(n.split, sum)
d.mean = sapply(dn.split, sum) / sapply(n.split, sum)
plotEmpir(results,log='xy',type='o')
lines(d.mean, e.mean, type='o', pch=19)
      

