## $Id:$
## Purpose: To generate the empirical SAD files

## read in community matrix files
shrtnames = c('bci','cocoli1','cocoli2','cross','sherman1','sherman2','sherman3',
              'serp')
comms = vector("list",length=length(shrtnames))
names(comms) = shrtnames
for (i in seq_along(comms))
  comms[[i]] = read.csv(paste0('./data/',shrtnames[i],'_comms.csv'))

## calculate empirical SAD
sadAbs = sapply(comms, function(x) 
                       apply(x[x$grain == unique(x$grain)[1], -(1:3)], 2, sum))

## export SAD files
for (i in seq_along(comms)) {
  write.table(matrix(sort(sadAbs[[i]], dec=TRUE), nrow=1), 
              file=paste0('./data/',shrtnames[i],'_sad.csv'),sep=',',
              row.names=FALSE,col.names=FALSE)
}