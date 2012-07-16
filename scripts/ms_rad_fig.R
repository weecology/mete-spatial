## manuscript and presentation style RAD figures

setwd('~/maxent/spat')

trop = c('bci', 'cocoli1', 'cocoli2', 'sherman1', 'sherman2', 'sherman3',
         'luquillo')
oakhick = c('oosting', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan', 
            'bigoak')
pine = c('graveyard', 'landsend')
savanna = 'cross'
evergreen = 'ferp'
grass = 'serp'

shrtnames = c(trop, oakhick, pine, savanna, evergreen, grass)
habindex = c(rep(1, length(trop)), rep(2, length(oakhick)), 
             rep(3, length(pine)), 4:6)
dat = vector("list", length=length(shrtnames))
names(dat) = shrtnames
for (i in seq_along(dat))
  dat[[i]] = read.csv(paste('./data/', shrtnames[i], '_sad.csv', sep=''),
                      header=FALSE)

shrtnm = read.table('./data/shrtnames.txt', colClasses='character')
N = read.table('./data/N_vals.txt')
S = read.table('./data/S_vals.txt')
N = as.numeric(N[match(names(dat), shrtnm)])
S = as.numeric(S[match(names(dat), shrtnm)])

freq = sapply(1:length(dat),  function(x) as.numeric(dat[[x]]) / N[x])
names(freq) = names(dat)
pdf('./figs/empirical_rads.pdf', width=7 * 2, height=7)
  col = terrain.colors(8)
  par(mfrow=c(1, 2))
  plot(1:S[1], as.numeric(freq[[1]]), ylim=range(freq), type='n', log='y', 
       xlab='rank', ylab='log Relative Freq.')
  for (i in seq_along(dat)) 
    lines(1:S[i], as.numeric(freq[[i]]), col=col[habindex[i]], lwd=2)
  legend('topright', c('tropical', 'oak-hickory', 'pine', 'oak savanna',
                       'mixed-evergreen', 'grassland'), col=col, lwd=2, bty='n')
  plot(1:S[1], as.numeric(freq[[1]]), ylim=range(freq), type='n', log='xy', 
       xlab='log rank', ylab='log Relative Freq.')
  for (i in seq_along(dat))
    lines(1:S[i], as.numeric(freq[[i]]), col=col[habindex[i]], lwd=2)
dev.off()


## for presentation
pdf('./figs/loglograd.pdf')
  col = terrain.colors(8)
  par(mfrow=c(1, 1))
  lwd=3
  plot(1:S[1], as.numeric(freq[[1]]), ylim=range(freq), type='n', log='xy', 
       xlab='', ylab='', frame.plot=F, axes=F)
  axis(side=1, cex.axis=1.5, lwd=3, padj=.25)
  axis(side=2, cex.axis=1.5, lwd=3)
  for (i in seq_along(dat))
    lines(1:S[i], as.numeric(freq[[i]]), col=col[habindex[i]], lwd=lwd)
  legend('topright', c('tropical', 'oak-hickory', 'old-field pine', 'oak savanna',
                       'mixed evergreen', 'grassland'), 
         col=col, lwd=lwd+2, cex=1.25, bty='n')
dev.off()
