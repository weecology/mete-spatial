## Author: Dan McGlinn
## Purpose: to generate a table that summarizes characteristics of the empirical
## datasets that are employed in the distance decay study

setwd('/home/danmcglinn/maxent/spat')
shrtnames = c('bci','cocoli1','cocoli2','sherman1','sherman2','sherman3','serp')
dat = vector("list",length=7)
names(dat) = shrtnames
for(i in seq_along(dat))
  dat[[i]] = read.csv(paste('./data/',shrtnames[i],'_sad.csv',sep=''),
                      header=FALSE)

datnames = c('BCI',rep('Cocoli',2),rep('Sherman',3),'Serpentine')
shape = c(rep('rectangle',5),rep('square',2))
area = c(50,rep(200*100*1e-3,4),140^2*1e-3,64*1e-3)
S = unlist(lapply(dat,length))
N = unlist(lapply(dat,sum))

datSummary = data.frame(datnames,shape,area,S,N)
write.csv(datSummary,file='./empir_data_summary.csv',row.names=FALSE)

freq = sapply(1:length(dat), function(x) as.numeric(dat[[x]]) / N[x])

pdf('empirical_rads.pdf')
par(mfrow=c(1,2))
plot(1:S[1],as.numeric(freq[[1]]),ylim=range(freq),type='n',log='y',
     xlab='rank',ylab='log Relative Freq.')
for(i in seq_along(dat))
  lines(1:S[i],as.numeric(freq[[i]]),col=i,lwd=2)
legend('topright',names(dat),col=1:length(dat),lwd=2,bty='n')
plot(1:S[1],as.numeric(freq[[1]]),ylim=range(freq),type='n',log='xy',
     xlab='log rank',ylab='log Relative Freq.')
for(i in seq_along(dat))
  lines(1:S[i],as.numeric(freq[[i]]),col=i,lwd=2)
dev.off()

