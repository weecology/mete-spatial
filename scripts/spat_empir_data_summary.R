## Author: Dan McGlinn
## Purpose: to generate a table that summarizes characteristics of the empirical
## datasets that are employed in the distance decay study

setwd('~/maxent/spat')
shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo')
dat = vector("list", length=length(shrtnames))
names(dat) = shrtnames
for (i in seq_along(dat))
  dat[[i]] = read.csv(paste0('./data/', shrtnames[i], '_sad.csv'), header=FALSE)

datnames = c('BCI', rep('Cocoli', 2), 'Crosstimbers', rep('Sherman', 3),
             'Serpentine', 'Oosting', 'FERP','Luquillo')
shape = c(rep('rectangle', 3), 'square', rep('rectangle', 2), rep('square', 3), 
          rep('rectangle', 2))
area = c(50, rep(200*100*1e-4, 2), 4, rep(200*100*1e-4, 2), 140^2*1e-4, 64*1e-4,
         160^2*1e-4, 150*300*1e-4, 250*500*1e-4)
S = unlist(lapply(dat, length))
N = unlist(lapply(dat, sum))

datSummary = data.frame(datnames, shape, area, S, N)
write.csv(datSummary, file='empir_data_summary.csv', row.names=FALSE)

freq = sapply(1:length(dat),  function(x) as.numeric(dat[[x]]) / N[x])

pdf('./figs/empirical_rads.pdf')
  par(mfrow=c(1, 2))
  plot(1:S[1], as.numeric(freq[[1]]), ylim=range(freq), type='n', log='y', 
       xlab='rank', ylab='log Relative Freq.')
  for (i in seq_along(dat))
    lines(1:S[i], as.numeric(freq[[i]]), col=i, lwd=2)
  legend('topright', names(dat), col=1:length(dat), lwd=2, bty='n')
  plot(1:S[1], as.numeric(freq[[1]]), ylim=range(freq), type='n', log='xy', 
       xlab='log rank', ylab='log Relative Freq.')
  for (i in seq_along(dat))
    lines(1:S[i], as.numeric(freq[[i]]), col=i, lwd=2)
dev.off()