## Author: Dan McGlinn
## Purpose: to generate a table that summarizes characteristics of the empirical
## datasets that are employed in the distance decay study

setwd('~/maxent/spat')
shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

dat = vector("list", length=length(shrtnames))
names(dat) = shrtnames
for (i in seq_along(dat))
  dat[[i]] = read.csv(paste('./data/', shrtnames[i], '_sad.csv', sep=''),
                      header=FALSE)

datnames = c('BCI', rep('Cocoli', 2), 'Crosstimbers', rep('Sherman', 3),
             'Serpentine', 'Oosting', 'FERP', 'Luquillo', 'Graveyard',
             'Landsend', 'Rocky', 'Bormann', 'Wood Bridge', 'Bald Mnt', 
             'Bryan', 'Big Oak')

shape = c(rep('rectangle', 3), 'square', rep('rectangle', 2), rep('square', 3), 
          rep('rectangle', 2), 'square', 'rectangle', rep('square', 3), 
          rep('rectangle', 3))

area = c(50 * 1e4, rep(200 * 100, 2), 4 * 1e4, rep(200 * 100, 2), 140^2, 64,
         160^2, 150 * 300, 250 * 500, 100^2, 130 * 65, 120^2, 140^2, 71^2,
         50 * 100, 185 * 92.5, 200 * 100) * 1e-4

S = unlist(lapply(dat, length))
N = unlist(lapply(dat, sum))

bisect = read.table('./data/bisect.txt', sep=' ', header=F)
bisect_fine = bisect[match(shrtnames, bisect[,1]), 2]
bisect_coarse = bisect[match(shrtnames, bisect[,1]), 3]

## output overall summary
datSummary = data.frame(datnames, shrtnames, shape, area_ha = area,
                        bisect_fine, bisect_coarse, 
                        grain_fine_m2 = area / 2^bisect_fine * 1e4,
                        grain_coarse_m2 = area / 2^bisect_coarse * 1e4,
                        S, N)
write.csv(datSummary, file='empir_data_summary.csv', row.names=FALSE)

## output shrtnames, S, N, and bisectoin info for the simulation routines
write.table(matrix(shrtnames, nrow=1), file=file.path('./data', 'shrtnames.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(matrix(S, nrow=1), file=file.path('./data', 'S_vals.txt'), sep=' ', 
            row.names=FALSE, col.names=FALSE)
write.table(matrix(N, nrow=1), file=file.path('./data', 'N_vals.txt'), sep=' ', 
            row.names=FALSE, col.names=FALSE)
write.table(matrix(bisect_fine, nrow=1), file=file.path('./data', 'bisect_fine.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE)
write.table(matrix(bisect_coarse, nrow=1), file=file.path('./data', 'bisect_coarse.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE)

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