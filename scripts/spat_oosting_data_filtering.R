## Purpose: to clean up the 19980 Oostings tree dataset for analysis of 
## 1990 data is used because its pre Hurrican Fran
## spatial biodiversity patterns. 
## Metadata: http://esapubs.org/archive/ecol/E088/162/metadata.htm

setwd('~/maxent/spat')

oost = read.delim('./data/raw_data/Oosting_Trees_1998.txt',
                  colClasses = 'character')

## replace the '.' symbols with NA's so that the data can be processed

oost = apply(oost, 2, function(x) ifelse(x == '.', NA, x))

dat = data.frame(ID = as.numeric(oost[,1]),
                 CODE = as.factor(oost[,2]),
                 NX = as.numeric(oost[,16]),
                 NY = as.numeric(oost[,17]),
                 D90 = as.numeric(oost[,5]),
                 C90 = as.numeric(oost[,6]))

## only keep the individual trees with a condition code of 1: alive and reasonably intact 
dat = subset(dat, C90 == 1)

write.csv(dat, file = './data/filtered_data/oosting_trees_1990_filtered.csv',
          row.names=F)

