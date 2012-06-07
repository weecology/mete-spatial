setwd('~/maxent/spat')

## read in comms file for the 1998 census of the crosstimbers
comms = read.csv('./data/cross_comms.csv', header=TRUE)
grains = unique(comms$grain)

## calculate empirical SAD/RAD
sadAbs = apply(comms[comms$grain == grains[1], -(1:3)], 2, sum)
write.table(matrix(sort(sadAbs, dec=TRUE), nrow=1), 
            file='./data/cross_sad.csv', sep=',', row.names=FALSE, 
            col.names=FALSE)

sadAvg = sapply(grains, function(x) apply(comms[comms$grain == x, -(1:3)], 2, mean))

## graphically explore SAD patterns
S= ncol(comms)-3
par(mfrow=c(1, 2))
plot(1:S, sort(sadAvg[ , 1], dec=TRUE), ylim=range(sadAvg), log='y', type='n')
for (i in seq_along(grains))
  lines(1:S, sort(sadAvg[ , i], dec=TRUE), type='o', col=i)
plot(1:S, sort(sadAvg[ , 1], dec=TRUE), ylim=range(sadAvg), log='xy', type='n')
for (i in seq_along(grains))
  lines(1:S, sort(sadAvg[ , i], dec=TRUE), type='o', col=i)


