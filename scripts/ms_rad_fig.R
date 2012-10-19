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

col = terrain.colors(8)
cols2 = c(104,128, 506, 139, 512, 92)
mcols2 = colors()[cols2]
col = mcols2


col = c('#b5f1b5', '#8bd3b1', '#7ee17e', '#bbeb80', '#66af66', '#749882', '#b6a092',
        '#facd84', '#faf591', '#d3e5ff', '#a9dff0', '#ebc2e2', '#d0cbf8', '#a1a9d3')
plot(1:length(col), col=col, pch=19, cex=3)
text(1:length(col), labels=1:length(col))

## trop, oak-hick, pine, savanna, evergreen, grass
indices = c(5, 11, 8, 7, 1, 12)
col = col[indices]
col =c("#66af66", "#a9dff0", "#facd84", "#b6a092", "#b5f1b5", "#ebc2e2")
col = c('forestgreen', 'medium purple', rev(colorschemes$Categorical.12[c(2, 4, 5, 8)]))
col = c("forestgreen", "#1AB2FF", "medium purple",
        "#B2FF8C", "#FFFF33", "#FF8000")
col = c("forestgreen", "#1AB2FF", "medium purple", 'saddlebrown', "#B2FF8C", 
        "#FF8000")
par(mfrow=c(1, 1))
lwd=3
plot(1:S[1], as.numeric(freq[[1]]), ylim=range(freq), type='n', log='xy', 
     xlab='', ylab='', frame.plot=F, axes=F)
axis(side=1, cex.axis=1.5, lwd=3, padj=.25)
axis(side=2, cex.axis=1.5, lwd=3)
for (i in seq_along(dat))
  lines(1:S[i], as.numeric(freq[[i]]), col=col[habindex[i]], lwd=lwd)

plot(1:10, 1:10, type='n', xlab='', ylab='', axes=F, frame.plot=F)
legend('center', c('tropical', 'oak-hickory', 'pine', 'oak savanna',
                     'mixed evergreen', 'grassland'), 
       col=col, lwd=lwd+4, cex=2, bty='n')


