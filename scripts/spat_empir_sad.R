## Purpose: To generate the empirical SAD files

print('Computing empirical SAD, ...')

## read in community matrix files
shrtnames = as.character(as.matrix(read.table('./data/shrtnames.txt')))

comms = vector("list", length=length(shrtnames))
names(comms) = shrtnames
for (i in seq_along(comms))
  comms[[i]] = read.csv(paste('./data/', shrtnames[i], '_comms.csv', sep=''))

## calculate empirical SAD
sadAbs = sapply(comms, function(x) 
                       apply(x[x$grain == unique(x$grain)[1], -(1:3)], 2, sum))

## drop species with zero abundances
sadAbs = sapply(sadAbs, function(x) x[x > 0])

## export SAD files
for (i in seq_along(comms)) {
  write.table(matrix(sort(sadAbs[[i]], dec=TRUE), nrow=1),  
              file=paste('./data/', shrtnames[i], '_sad.csv', sep=''), sep=',', 
              row.names=FALSE, col.names=FALSE)
}

print('Computing empirical SAD, complete!')
