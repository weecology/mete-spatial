## Purpose: to graphically summarize the simulated DDR results for each
## empirical dataset
setwd('~/maxent/spat')
source('./scripts/spat_sim_vario_func.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

files = get_file_names(path = './sorensen/', 
                      strings = list(shrtnames, 'C200', 'sherman3'), 
                      invert = c(F, F, T))

longnames = sub('sorensen_', '', files)
longnames = sub('_abu.Rdata', '', longnames)
longnames = unique(sub('_binary.Rdata', '', longnames))

simSorBin = getResults(longnames, 'sorensen', 'binary', sim_result=TRUE)
simSorAbu = getResults(longnames, 'sorensen', 'abu', sim_result=TRUE)

simSorBinAvg = avgSimResults(simSorBin, 'sorensen')
simSorAbuAvg = avgSimResults(simSorAbu, 'sorensen', null=FALSE)

simVarBin = getResults(longnames, 'varWithin', 'binary', sim_result=TRUE)
simVarAbu = getResults(longnames, 'varWithin', 'abu', sim_result=TRUE)

simVarBinAvg = avgSimResults(simVarBin, 'varWithin')
simVarAbuAvg = avgSimResults(simVarAbu, 'varWithin')

## change names to better match empirical names
for (i in seq_along(longnames)) {
  newnames = c('Comm', 'Dist', 'N', 'Med.lo', 'Med', 'Med.hi', 'Avg.lo', 'Avg',
               'Avg.hi', 'Exp.lo', 'Exp', 'Exp.hi')
  names(simSorBinAvg[[i]]) = newnames
  names(simSorAbuAvg[[i]]) = newnames[-(10:12)]
  names(simVarBinAvg[[i]]) = newnames
  names(simVarAbuAvg[[i]]) = newnames
}

## export as files
#save(simSorBinAvg, simSorAbuAvg, file = './sorensen/simEmpirSorAvg.Rdata')
#save(simVarBinAvg, simVarAbuAvg, file = './varWithin/simEmpirVarAvg.Rdata')
load('./sorensen/simEmpirSorAvg.Rdata')
load('./varWithin/simEmpirVarAvg.Rdata')

## split results that were generated from log-series and those that were fixed abu
logser = grep('empirSAD', longnames, invert=TRUE)
fixed = grep('empirSAD', longnames, invert=FALSE)
simSorBinLogSer = simSorBinAvg[logser]
simSorBinFixed = simSorBinAvg[fixed]

simSorAbuLogSer = simSorAbuAvg[logser]
simSorAbuFixed = simSorAbuAvg[fixed]

simVarBinLogSer = simVarBinAvg[logser]
simVarBinFixed = simVarBinAvg[fixed]

simVarAbuLogSer = simVarAbuAvg[logser]
simVarAbuFixed = simVarAbuAvg[fixed]

## merge cocoli 1 & 2 and sherman 1 & 2
simSorBinLogSer = merge_drop(simSorBinLogSer)
simSorBinFixed = merge_drop(simSorBinFixed )

simSorAbuLogSer = merge_drop(simSorAbuLogSer)
simSorAbuFixed = merge_drop(simSorAbuFixed)

simVarBinLogSer = merge_drop(simVarBinLogSer)
simVarBinFixed = merge_drop(simVarBinFixed)

simVarAbuLogSer = merge_drop(simVarAbuLogSer)
simVarAbuFixed = merge_drop(simVarAbuFixed)

## update shrnames
shrtnames = shrtnames[-match(c('cocoli2','sherman2', 'sherman3'), shrtnames)] 

## export results
save(simSorAbuLogSer, simSorAbuFixed, simSorBinLogSer, simSorBinFixed,
     simVarAbuLogSer, simVarAbuFixed, simVarBinLogSer, simVarBinFixed,
     file = './simulated_empirical_results.Rdata')

## load the Results
load('./simulated_empirical_results.Rdata')


## examine a single simulated distance decay pattern----------------------------
tmp = simSorAbuLogSer[3]
col = colorRampPalette(c('dodgerblue', 'red'))(5)
range(tmp[[1]]$Avg)
grains = unique(tmp[[1]]$Comm)

plot(1:10, type='n', xlab='', ylab='', xlim=c(0, 110), 
     ylim = c(0, 0.4) , frame.plot=F, axes=F, log='')
axis(side=1, cex.axis=1.25, lwd=2)
axis(side=2, cex.axis=1.25, lwd=2)
plotEmpir(tmp, 'average', log='xy', title=F, 
          quants=F, col=col, lwd=5, add=TRUE, type='p', pch=19)
for (g in seq_along(grains)) {
  mod = lm(log(Avg) ~ log(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  lines(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]], exp(predict(mod)),
        col=col[g], lwd=4)  
}

###
plot(1:10, type='n', xlab='', ylab='', xlim=c(1.5, 120), 
     ylim = c(0.006, 0.40) , frame.plot=F, axes=F, log='xy')
axis(side=1, cex.axis=1.75, padj=.5, lwd=8)
axis(side=2, cex.axis=1.75, lwd=8, at=c(0.01, 0.02, 0.05, 0.10, 0.20, 0.40))
plotEmpir(tmp, 'average', log='xy', title=F, 
          quants=F, colcol, lwd=5, add=TRUE, type='p', cex=1.5, pch=19)
for (g in seq_along(grains)) {
  mod = lm(log(Avg) ~ log(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  lines(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]], exp(predict(mod)),
        col=col[g], lwd=6)  
}
##
plot(log2(Avg) ~ log2(Dist), data=tmp[[1]], type='n', xlab='', ylab='',
     frame.plot=F, axes=F)
axis(side=1, cex.axis=1.25, lwd=2)
axis(side=2, cex.axis=1.25, lwd=2)
for (g in seq_along(grains)) {
  mod = lm(log2(Avg) ~ log2(Dist), data=tmp[[1]], subset=Comm == grains[g],
           weights = N)
  points(log2(Avg) ~ log2(Dist), data=tmp[[1]], subset=Comm == grains[g],
         col=col[g], pch=19, lwd=5)
  lines(log2(tmp[[1]]$Dist[tmp[[1]]$Comm == grains[g]]), predict(mod),
        col=col[g], lwd=4)  
}


## export pdf summaries of the empirical and simulated results------------------
mk_summary = function() {
    par(mfrow=c(2,3))
    comp_results(site,
                 list(results1),
                 list("quants=F, expected=TRUE, sub='RandomPlacement',
                       metric='average', ylim=ylim, col=col, log=log_tran"))
    comp_results(site,
                 list(results1, results2),
                 list("quants=F, sub='Fixed'",
                      "quants=T, lwd=2, lty=2"),
                 metric='average', ylim=ylim, col=col, log=log_tran)
    comp_results(site,
                 list(results1, results3),
                 list("quants=F, sub='LogSer'",
                      "quants=T, lwd=2, lty=2"),
                 metric='average', ylim=ylim,col=col, log=log_tran)
    ###
    comp_results(site,
                 list(results1),
                 list("quants=F, ylim=ylim,
                       col=col, log=log_tran"))
    comp_results(site,
                 list(results1, results2),
                 list("quants=T, sub='AbuFixed'",
                      "quants=T, lwd=2, lty=2"),
                 metric='median', ylim=ylim, col=col, log=log_tran)
    comp_results(site,
                 list(results1, results3),
                 list("quants=T, sub='LogSer'",
                      "quants=T, lwd=2, lty=2"),
                 metric='median', ylim=ylim,col=col, log=log_tran)
}

##------------------------------------------------------------------------------
results1 = empirSorAbu
results2 = simSorAbuFixed
results3 = simSorAbuLogSer
metric = 'sorensen'
dataType = 'abu'
log_tran='xy'
ylim=list(c(0.1,1))
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------

results1 = empirSorBin
results2 = simSorBinFixed
results3 = simSorBinLogSer
metric = 'sorensen'
dataType = 'bin'
log_tran='xy'
ylim=list(c(0.1,1))
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------
results1 = empirVarAbu
results2 = simVarAbuFixed
results3 = simVarAbuLogSer
metric = 'varWithin'
dataType = 'abu'
log_tran='xy'
ylim=NULL
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
##------------------------------------------------------------------------------

results1 = empirVarBin
results2 = simVarBinFixed
results3 = simVarBinLogSer
metric = 'varWithin'
dataType = 'bin'
log_tran='xy'
ylim=NULL
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
  for (i in seq_along(results2)) {
    site = names(results1)[i]
    col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                                '$Comm)) + 2)', sep='')))
    mk_summary()
  }
dev.off()
