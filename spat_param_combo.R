
##what parameter combinations to consider


S <- round(10^seq(log10(10),log10(100),length.out=20))

N <- round(10^seq(log10(120),log10(5e5),length.out=20))


grid <- expand.grid(N,S)
par(mfrow=c(1,2))
plot(grid)
plot(grid[,2]/grid[,1])
ratios<-round(sort(unique(grid[,2]/grid[,1])),3)

grid.uni <- grid[unique(match(ratios,round(grid[,2]/grid[,1],3))),]
dim(grid.uni)

par(mfrow=c(2,2))
plot(grid,log='xy',xlab='log N',ylab='log S')
plot(sort(grid[,2]/grid[,1],decreasing = TRUE), xlab='rank',ylab='S/N')
plot(grid.uni,log='xy',xlab='log N',ylab='log S')
plot(sort(grid.uni[,2]/grid.uni[,1],decreasing = TRUE), xlab='rank',ylab='S/N')


