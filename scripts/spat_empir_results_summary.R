## Purpose: to graphically summarize the empirical results
setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

empirBin = getResults(shrtnames,'sorensen','binary')
empirAbu = getResults(shrtnames,'sorensen','abu')
empirSorBin = reshapeResults(empirBin,'sorensen')
empirSorAbu = reshapeResults(empirAbu,'sorensen')
## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirSorBin = merge_drop(empirSorBin)
empirSorAbu = merge_drop(empirSorAbu)

empirBin = getResults(shrtnames,'varWithin','binary')
empirAbu = getResults(shrtnames,'varWithin','abu')
empirVarBin = reshapeResults(empirBin,'varWithin')
empirVarAbu = reshapeResults(empirAbu,'varWithin')
## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirVarBin = merge_drop(empirVarBin)
empirVarAbu = merge_drop(empirVarAbu)

## export results to file
save(empirSorBin, file='./sorensen/empirSorBin.Rdata')
save(empirSorAbu, file='./sorensen/empirSorAbu.Rdata')
save(empirVarBin, file='./varWithin/empirVarBin.Rdata')
save(empirVarAbu, file='./varWithin/empirVarAbu.Rdata')

## examine results
par(mfrow=c(4,4))
plotEmpir(empirSorAbu, type='o')
plotEmpir(empirSorAbu, log='xy', type='o')
plotEmpir(empirSorAbu,log='xy', type='o',quants=T)

par(mfrow=c(4,4))
plotEmpir(empirVarAbu,log='x', type='o')
plotEmpir(empirVarAbu,log='xy', type='o',quants=T)

##
commName = 'cocoli1'
par(mfrow=c(1,2))
plotEmpir(empirSorAbu[commName],log='xy', type='o', quants=T)
plotEmpir(empirVarAbu[commName],log='xy', type='o', quants=F)

## mean vs median
plotEmpir(empirSorAbu[1], metric='median', log='xy')
plotEmpir(empirSorAbu[1], metric='average', log='xy', lty=2, add=T)

plotEmpir(empirSorBin[1], metric='median', log='xy')
plotEmpir(empirSorBin[1], metric='average', log='xy', lty=2, add=T)

plotEmpir(empirVarAbu[1],  metric='median', log='xy')
plotEmpir(empirVarAbu[1], metric='average', log='xy', lty=2, add=T)

plotEmpir(empirVarBin[1],  metric='median', log='xy')
plotEmpir(empirVarBin[1], metric='average', log='xy', lty=2, add=T)

## generate pdfs
for (metric in c('Sor', 'Var')) {
  for (type in c('Bin', 'Abu')) {
    for (log_it in c(FALSE, TRUE)) {
      if (log_it)
        prefix = './figs/spat_empir_loglog_'
      else
        prefix = './figs/spat_empir_'
      pdf(paste(prefix, metric, '_', type, '_curves.pdf', sep=''),
          width = 7 * 2, height = 7 * 2)
        par(mfrow=c(4,4))
        data = eval(parse(text=paste('empir', metric, type, sep='')))
        if (log_it) {
          plotEmpir(data, quants=FALSE, type='o', log='xy')
          plotEmpir(data, quants=TRUE, type='o', log='xy')
        }  
        else {
          plotEmpir(data, quants=FALSE, type='o')
          plotEmpir(data, quants=TRUE, type='o')
        }
        rm(data)
      dev.off()
    }  
  }  
}
