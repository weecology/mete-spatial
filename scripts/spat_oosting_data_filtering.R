## Purpose: to clean up the 1998 Oostings tree dataset for analysis of 
## 1990 data is used because its pre Hurrican Fran
## spatial biodiversity patterns. 
## Metadata: http://esapubs.org/archive/ecol/E088/162/metadata.htm

setwd('~/maxent/spat')

dat = read.delim('./data/raw_data/Oosting_Trees_1998.txt',
                  colClasses = 'character')

## replace the '.' symbols with NA's so that the data can be processed

dat = apply(dat, 2, function(x) ifelse(x == '.', NA, x))

oost = data.frame(ID = as.numeric(dat[,1]),
                  CODE = as.factor(dat[,2]),
                  NX = as.numeric(dat[,16]),
                  NY = as.numeric(dat[,17]),
                  D90 = as.numeric(dat[,5]),
                  C90 = as.numeric(dat[,6]))

## only keep the individual trees with a condition code of 1: alive and reasonably intact 
oost = subset(oost, C90 == 1)

write.csv(oost, file = './data/filtered_data/oosting_trees_1990_filtered.csv',
          row.names=F)

