setwd('~/maxent/spat')

## read in comms file for the 1999 census (i.e.,  census 3) of sherman
comms1 = read.csv('./data/sherman1_comms.csv')
comms2 = read.csv('./data/sherman2_comms.csv')
comms3 = read.csv('./data/sherman3_comms.csv')

## there are two 200 x 100 m quadrats and one 140 x 140 m quad
grains12 = unique(comms1$grain)
grains3 = unique(comms3$grain)

## calculate empirical SAD/RAD for each
sadAbs1 = apply(comms1[comms1$grain == grains12[1], -(1:3)], 2, sum)
sadAbs2 = apply(comms2[comms2$grain == grains12[1], -(1:3)], 2, sum)
sadAbs3 = apply(comms3[comms3$grain == grains3[1], -(1:3)], 2, sum)

## drop species in each with abundances of zero
sadAbs1 = sadAbs1[sadAbs1 > 0]
sadAbs2 = sadAbs2[sadAbs2 > 0]
sadAbs3 = sadAbs3[sadAbs3 > 0]

write.table(matrix(sort(sadAbs1, dec=TRUE), nrow=1), 
            file='./data/sherman1_sad.csv', sep=',', row.names=FALSE, 
            col.names=FALSE)
write.table(matrix(sort(sadAbs2, dec=TRUE), nrow=1), 
            file='./data/sherman2_sad.csv', sep=',', row.names=FALSE, 
            col.names=FALSE)
write.table(matrix(sort(sadAbs3, dec=TRUE), nrow=1), 
            file='./data/sherman3_sad.csv', sep=',', row.names=FALSE, 
            col.names=FALSE)
            
sadAvg = sapply(grains, function(x) apply(comms[comms$grain == x, -(1:3)], 2, mean))

## graphically explore SAD patterns
S= ncol(comms) - 3
par(mfrow=c(1, 2))
plot(1:S, sort(sadAvg[ , 1] + 1, dec=TRUE), ylim=range(sadAvg) + 1, log='y',
     type='n')
for (i in seq_along(grains))
  lines(1:S, sort(sadAvg[ , i] + 1, dec=TRUE), type='o', col=i)
plot(1:S, sort(sadAvg[ , 1] + 1, dec=TRUE), ylim=range(sadAvg) + 1, log='xy',
     type='n')
for (i in seq_along(grains))
  lines(1:S, sort(sadAvg[ , i] + 1, dec=TRUE), type='o', col=i)

