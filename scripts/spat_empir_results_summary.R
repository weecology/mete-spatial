## Purpose: to graphically summarize the empirical results
source('./scripts/spat_functions.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

empirBin = getResults(shrtnames, 'sorensen', 'binary', bisect=TRUE)
empirAbu = getResults(shrtnames, 'sorensen', 'abu', bisect=TRUE)
empirSorBin = reshapeResults(empirBin, 'sorensen', bisect=TRUE)
empirSorAbu = reshapeResults(empirAbu, 'sorensen', perm_null=TRUE, bisect=TRUE)

## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirSorBin = merge_drop(empirSorBin)
empirSorAbu = merge_drop(empirSorAbu)

## export results to file
save(empirSorBin, file='./sorensen/empirSorBin_bisect.Rdata')
save(empirSorAbu, file='./sorensen/empirSorAbu_bisect.Rdata')

## examine results
par(mfrow=c(2,4))
plotEmpir(empirSorAbu, type='o')
plotEmpir(empirSorAbu, log='xy', type='o')
plotEmpir(empirSorAbu,log='xy', type='o',quants=T)

##
commName = 'serp'
par(mfrow=c(1,1))
plotEmpir(empirSorAbu[commName], metric='average', type='o')
plotEmpir(empirSorAbu[commName],log='xy', type='o', quants=T)

## mean vs median
plotEmpir(empirSorAbu[1], metric='median', log='xy')
plotEmpir(empirSorAbu[1], metric='average', log='xy', lty=2, add=T)

plotEmpir(empirSorBin[1], metric='median', log='xy')
plotEmpir(empirSorBin[1], metric='average', log='xy', lty=2, add=T)

## generate pdfs
for (metric in c('Sor')) {
  for (type in c('Bin', 'Abu')) {
    for (log_it in c(FALSE, TRUE)) {
      if (log_it)
        prefix = './figs/spat_empir_loglog_bisect_'
      else
        prefix = './figs/spat_empir_bisect_'
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
