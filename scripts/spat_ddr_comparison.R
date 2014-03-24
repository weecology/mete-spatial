## Purpose: to compare the fit of the simulated ddr to the empirical ddr

source('./scripts/spat_functions.R')

## load the simulated and empirical DDR results---------------------------------
load('./simulated_empirical_results_bisect.Rdata')
load('./sorensen/empirSorBin_bisect.Rdata') 
load('./sorensen/empirSorAbu_bisect.Rdata') 

## compute residuals between observed and expected-------------------------------
resSorAbuFixed = get_ddr_resid(empirSorAbu, simSorAbuFixed)
resSorAbuLogSer = get_ddr_resid(empirSorAbu, simSorAbuLogSer)

resSorBinFixed = get_ddr_resid(empirSorBin, simSorBinFixed)
resSorBinLogSer = get_ddr_resid(empirSorBin, simSorBinLogSer)


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
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_bisect_comp.pdf', sep=''),
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
pdf(paste('./figs/empir_sim_', metric, '_', dataType, '_bisect_comp.pdf', sep=''),
    width=7*2, height=7*1.5)
for (i in seq_along(results2)) {
  site = names(results1)[i]
  col = eval(parse(text=paste('terrain.colors(length(unique(results1$', site,
                              '$Comm)) + 2)', sep='')))
  mk_summary()
}
dev.off()


