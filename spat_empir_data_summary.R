## Author: Dan McGlinn
## Purpose: to generate a table that summarizes characteristics of the empirical
## datasets that are employed in the distance decay study

setwd('/home/danmcglinn/maxent/spat')

dat = list()
dat[[1]] = read.csv('./data/bci_sad.csv',header=FALSE)
dat[[2]] = read.csv('./data/cocoli_sad_1.csv',header=FALSE)
dat[[3]] = read.csv('./data/cocoli_sad_2.csv',header=FALSE)
dat[[4]] = read.csv('./data/sherman_sad_1.csv',header=FALSE)
dat[[5]] = read.csv('./data/sherman_sad_2.csv',header=FALSE)
dat[[6]] = read.csv('./data/sherman_sad_3.csv',header=FALSE)
dat[[7]] = read.csv('./data/serpentine_sad.csv',header=FALSE)
names(dat) = c('bci','cocoli1','cocoli2','sherman1','sherman2','sherman3','serp')

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

