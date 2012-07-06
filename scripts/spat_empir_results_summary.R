## Purpose: to graphically summarize the empirical results
setwd('./maxent/spat')
source('./scripts/spat_sim_vario_func.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

merge_drop = function(results) {
  result_names = names(results)
  to_drop = match(c('cocoli2','sherman2','sherman3'), result_names)
  to_avg = c('Metric.avg', 'Metric.25', 'Metric.50', 'Metric.75', 'Exp', 'N')
  cocoli_results = (subset(results$cocoli1, select= -Comm) + 
                    subset(results$cocoli2, select= -Comm)) / 2
  cocoli_results = data.frame(cocoli_results, Comm=results$cocoli1[ , 'Comm'])
  sherman_results = (subset(results$sherman1, select= -Comm) + 
                     subset(results$sherman2, select= -Comm)) / 2
  sherman_results = data.frame(sherman_results, Comm=results$sherman1[ , 'Comm'])
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

empirBin = getResults(shrtnames,'varWithin','binary')
empirAbu = getResults(shrtnames,'varWithin','abu')
empirVarBin = reshapeResults(empirBin,'varWithin')
empirVarAbu = reshapeResults(empirAbu,'varWithin')
## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirVarBin = merge_drop(empirVarBin)
empirVarAbu = merge_drop(empirVarAbu)

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
