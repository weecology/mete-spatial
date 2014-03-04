## Purpose: to graphically summarize the simulated DDR results for each
## empirical dataset

print('Aggregating and exporting the simulated DDR results, ...')

source('./scripts/spat_functions.R')

shrtnames = as.character(read.table('./data/shrtnames.txt', colClasses='character')[1,])

files = get_file_names(path = './sorensen/', 
                      strings = list(shrtnames, 'C200', 'bisect', 'sherman3'), 
                      invert = c(F, F, F, T))

longnames = sub('sorensen_', '', files)
longnames = sub('_bisect_abu.Rdata', '', longnames)
longnames = unique(sub('_bisect_binary.Rdata', '', longnames))

simSorBin = getResults(longnames, 'sorensen', 'binary', bisect=TRUE, sim_result=TRUE)
simSorAbu = getResults(longnames, 'sorensen', 'abu', bisect=TRUE, sim_result=TRUE)

simSorBinAvg = avgSimResults(simSorBin, 'sorensen', bisect=TRUE)
simSorAbuAvg = avgSimResults(simSorAbu, 'sorensen', null=FALSE, bisect=TRUE)

## change names to better match empirical names
for (i in seq_along(longnames)) {
  newnames = c('Comm', 'Dist', 'N', 'Med.lo', 'Med', 'Med.hi', 'Avg.lo', 'Avg',
               'Avg.hi', 'Exp.lo', 'Exp', 'Exp.hi')
  names(simSorBinAvg[[i]]) = newnames
  names(simSorAbuAvg[[i]]) = newnames[-(10:12)]
}

## export the aggregrated raw results files
save(simSorBinAvg, simSorAbuAvg, file = './sorensen/simEmpirSorAvg_bisect.Rdata')

## split results that were generated from log-series and those that were fixed abu
logser = grep('empirSAD', longnames, invert=TRUE)
fixed = grep('empirSAD', longnames, invert=FALSE)
simSorBinLogSer = simSorBinAvg[logser]
simSorBinFixed = simSorBinAvg[fixed]

simSorAbuLogSer = simSorAbuAvg[logser]
simSorAbuFixed = simSorAbuAvg[fixed]

## merge cocoli 1 & 2 and sherman 1 & 2
simSorBinLogSer = merge_drop(simSorBinLogSer)
simSorBinFixed = merge_drop(simSorBinFixed )

simSorAbuLogSer = merge_drop(simSorAbuLogSer)
simSorAbuFixed = merge_drop(simSorAbuFixed)

## export the agregated merged results
save(simSorAbuLogSer, simSorAbuFixed, simSorBinLogSer, simSorBinFixed,
     file = './simulated_empirical_results_bisect.Rdata')

print('Aggregating and exporting the simulated DDR results, complete!')

