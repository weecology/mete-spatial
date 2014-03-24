## Purpose: to graphically summarize the empirical results

print('Aggregating and exporting the empirical DDR results, ...')

source('./scripts/spat_functions.R')

shrtnames = as.character(read.table('./data/shrtnames.txt', colClasses='character')[1,])

empirBin = getResults(shrtnames, 'sorensen', 'binary', bisect=TRUE)
empirAbu = getResults(shrtnames, 'sorensen', 'abu', bisect=TRUE, swap='indiv')

empirSorBin = reshapeResults(empirBin, 'sorensen', bisect=TRUE)
empirSorAbu = reshapeResults(empirAbu, 'sorensen', perm_null=TRUE, bisect=TRUE)

## Average cocoli1 & 2 and sherman 1 & 2 and drop sherman3
empirSorBin = merge_drop(empirSorBin)
empirSorAbu = merge_drop(empirSorAbu)

## export results to file
save(empirSorBin, file='./sorensen/empirSorBin_bisect.Rdata')
save(empirSorAbu, file='./sorensen/empirSorAbu_bisect.Rdata')

print('Aggregating and exporting the empirical DDR results, complete!')
