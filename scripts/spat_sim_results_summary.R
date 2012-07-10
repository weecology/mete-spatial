## Purpose: to graphically summarize the simulated simical results
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

simSorBin = getResults(longnames, 'sorensen', 'binary', TRUE)
simSorAbu = getResults(longnames, 'sorensen', 'abu', TRUE)

simSorBinAvg = avgSimResults(simSorBin, 'sorensen')
simSorAbuAvg = avgSimResults(simSorAbu, 'sorensen')

simVarBin = getResults(longnames, 'varWithin', 'binary', TRUE)
simVarAbu = getResults(longnames, 'varWithin', 'abu', TRUE)

simVarBinAvg = avgSimResults(simVarBin, 'varWithin')
simVarAbuAvg = avgSimResults(simVarAbu, 'varWithin')

## change names to better match empirical names
for (i in seq_along(longnames)) {
  newnames = c('Comm', 'Dist', 'N', 'Med.lo', 'Med', 'Med.hi', 'Avg.lo', 'Avg',
               'Avg.hi', 'Exp.lo', 'Exp', 'Exp.hi')
  names(simSorBinAvg[[i]]) = newnames
  names(simSorAbuAvg[[i]]) = newnames
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

## Now bring in empirical Results
load('./sorensen/empirSorBin.Rdata') 
load('./sorensen/empirSorAbu.Rdata') 
load('./varWithin/empirVarBin.Rdata')
load('./varWithin/empirVarAbu.Rdata')

## export pdf summaries of the empirical and simulated results
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

results1 = empirBinAbu
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
##------------------------------------------------------------------------------

## compute summary statics for empricial and simualted results
empirStatsSorAbuMed = getStats(empirSorAbu, 'median')
empirStatsSorAbuAvg = getStats(empirSorAbu, 'average')

empirStatsSorAbuMed$ferp = empirStatsSorAbuMed$ferp[,,,-6]
empirStatsSorAbuAvg$ferp = empirStatsSorAbuAvg$ferp[,,,-6]
  
simStatsSorAbuFixedMed = getStats(simSorAbuFixed, 'median')
simStatsSorAbuFixedAvg = getStats(simSorAbuFixed, 'average')

simStatsSorAbuLogSerMed = getStats(simSorAbuLogSer, 'median')
simStatsSorAbuLogSerAvg = getStats(simSorAbuLogSer, 'average')

##
empirStatsVarAbuMed = getStats(empirVarAbu, 'median')
empirStatsVarAbuAvg = getStats(empirVarAbu, 'average')

empirStatsVarAbuMed$ferp = empirStatsVarAbuMed$ferp[,,,-6]
empirStatsVarAbuAvg$ferp = empirStatsVarAbuAvg$ferp[,,,-6]
  
simStatsVarAbuFixedMed = getStats(simVarAbuFixed, 'median')
simStatsVarAbuFixedAvg = getStats(simVarAbuFixed, 'average')

simStatsVarAbuLogSerMed = getStats(simVarAbuLogSer, 'median')
simStatsVarAbuLogSerAvg = getStats(simVarAbuLogSer, 'average')

##-----------------------------------------------------------------------------

empirStatsSorBinMed = getStats(empirSorBin, 'median')
empirStatsSorBinAvg = getStats(empirSorBin, 'average')

empirStatsSorBinMed$ferp = empirStatsSorBinMed$ferp[,,,-6]
empirStatsSorBinAvg$ferp = empirStatsSorBinAvg$ferp[,,,-6]
  
simStatsSorBinFixedMed = getStats(simSorBinFixed, 'median')
simStatsSorBinFixedAvg = getStats(simSorBinFixed, 'average')

simStatsSorBinLogSerMed = getStats(simSorBinLogSer, 'median')
simStatsSorBinLogSerAvg = getStats(simSorBinLogSer, 'average')

##
empirStatsVarBinMed = getStats(empirVarBin, 'median')
empirStatsVarBinAvg = getStats(empirVarBin, 'average')

empirStatsVarBinMed$ferp = empirStatsVarBinMed$ferp[,,,-6]
empirStatsVarBinAvg$ferp = empirStatsVarBinAvg$ferp[,,,-6]
  
simStatsVarBinFixedMed = getStats(simVarBinFixed, 'median')
simStatsVarBinFixedAvg = getStats(simVarBinFixed, 'average')

simStatsVarBinLogSerMed = getStats(simVarBinLogSer, 'median')
simStatsVarBinLogSerAvg = getStats(simVarBinLogSer, 'average')


## generate pdfs of r2 values, exponential or pwr function better
pl_r2 = function() {
  for( i in 1:3) {
    r2pwr = unlist(sapply(stats[[i]], function(x) x['pwr', 'r2', 'wtr',]))
    r2exp = unlist(sapply(stats[[i]], function(x) x['exp', 'r2', 'wtr',]))
    brks = seq(0, 1, length.out=21)
    par(mfrow=c(1,2))
    hist(r2exp, xlim=c(0,1), breaks=brks, col='red')
    hist(r2pwr, xlim=c(0,1), breaks=brks, col='blue', add=TRUE)
    plot(r2pwr, r2exp)
    abline(a=0, b=1)
  }    
}

stats = list(empirStatsSorAbuAvg, simStatsSorAbuFixedAvg, simStatsSorAbuLogSerAvg)
pdf('./figs/exp_vs_pwr_avg_sor_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsSorAbuMed, simStatsSorAbuFixedMed, simStatsSorAbuLogSerMed)
pdf('./figs/exp_vs_pwr_med_sor_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsVarAbuAvg, simStatsVarAbuFixedAvg, simStatsVarAbuLogSerAvg)
pdf('./figs/exp_vs_pwr_avg_var_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  
##----------------------------------------------------------------------------
stats = list(empirStatsVarAbuMed, simStatsVarAbuFixedMed, simStatsVarAbuLogSerMed)
pdf('./figs/exp_vs_pwr_med_var_abu_empirical.pdf', width= 7 * 2, height=7) 
  pl_r2()
dev.off()  

## plot summary statistics

pl_coef = function(stats1, stats2, xlim=NULL, ylim=NULL) {
  if (is.null(xlim))
    xlim = range(unlist(sapply(stats1, function(x) range(x[mod, cof, meth, ], na.rm=T))))
  if (is.null(ylim))
    ylim = range(unlist(sapply(stats2, function(x) range(x[mod, cof, meth, ], na.rm=T))))
  plot(stats1$bci[mod, cof, meth, ], stats2$bci[mod, cof, meth, ], 
       type='n', xlim = xlim, ylim= ylim, main=paste(mod, meth, sep=', '),
       xlab=paste('Empir', cof, sep=' '), ylab=paste('METE', cof, sep=' '))
  for (i in seq_along(stats1)) {
    site = names(stats1)[i]
    index2 = which(grepl(site, names(stats2)))
    index3 = which(grepl(site, names(stats3)))
    if (length(index2) > 0)
      points(stats1[[i]][mod, cof, meth, ], 
             stats2[[index2]][mod, cof, meth, ], col='blue')
    if (length(index3) > 0)
      points(stats1[[i]][mod, cof, meth, ], 
             stats3[[index3]][mod, cof, meth, ], col='red')
  }
  abline(a=0, b=1)
}

pdf('./figs/coef_one_to_one_sor_abu_bin_avg_wtr.pdf', width=7 * 2, height=7*1.5)
  stats1 = empirStatsSorAbuAvg
  stats2 = simStatsSorAbuFixedAvg
  stats3 = simStatsSorAbuLogSerAvg
  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
##------------------------------
  stats1 = empirStatsSorBinAvg
  stats2 = simStatsSorBinFixedAvg
  stats3 = simStatsSorBinLogSerAvg
  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
dev.off()

##-----------------------------------------------------------------------------
pdf('./figs/coef_one_to_one_var_abu_bin_avg_wtr.pdf', width=7 * 2, height=7*1.5)
  stats1 = empirStatsVarAbuAvg 
  stats2 = simStatsVarAbuFixedAvg
  stats3 = simStatsVarAbuLogSerAvg
  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
##------------------------------
  stats1 = empirStatsVarBinAvg 
  stats2 = simStatsVarBinFixedAvg
  stats3 = simStatsVarBinLogSerAvg
  par(mfrow=c(2,2))
  mod = 'pwr'
  meth = 'wtr'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
  ##
  mod = 'exp'
  cof = 'b0'
  pl_coef(stats1, stats2)
  #
  cof = 'b1'
  pl_coef(stats1, stats2)
dev.off()  

################ Abu Sorensen
pdf('empir_sim_sor_arith.pdf')
par(mfrow=c(4,4))
ylims = rbind(c(0.001,.4),c(.001,.05),c(.001,.05),c(.001,.4),
              c(.001,.7),c(.001,.05),c(.001,.05),c(.001,.05))
for (i in seq_along(shrtnames)) {
  results1 = eval(parse(text=paste('empirSorAbu$', shrtnames[i], sep='')))           
  results2 = eval(parse(text=paste('simSorAbuLogSer$', shrtnames[i], sep=''))) 
  results3 = eval(parse(text=paste('simSorAbuFixed$', shrtnames[i], sep=''))) 
  unigrains = unique(results1$Comm)
  plot(Metric.50 ~ Dist, data = results1, ylim=c(0, 1), type='n', main=shrtnames[i])
  for (j in seq_along(unigrains)) {
    lines(Metric.50 ~ Dist, data = results1, subset= Comm == unigrains[j], 
         col=j, lty=1, lwd=2, type='o')
    lines(Med ~ Dist, data = results3, subset= Comm == unigrains[j],
          col=j, lty=2, lwd=2, type='l')
    if (shrtnames[i] != 'cross')
      lines(Med ~ Dist, data = results2, subset= Comm == unigrains[j],
            col=j, lty=3, lwd=2, type='l')
  }
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
ylims = rbind(c(0.02,1),c(.001,1),c(.001,1),c(.02,.4),
              c(.1,1),c(.003,1),c(.003,1),c(.003,1))
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
ylims = rbind(c(0.05,1),c(.001,1),c(.001,1),c(.03,.8),
              c(.2,1),c(.003,1),c(.003,1),c(.003,1))
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
