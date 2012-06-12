## Purpose: to clean up the 2007 Ferp woodly plant dataset for analysis of 
## spatial biodiversity patterns. 
## Metadata: http://ferp.ucsc.edu/ferp/about/using/

setwd('~/maxent/spat')

dat= read.delim('~/datasets/CTFSplots/FERP/PlotsDataReport.txt')

## drop dead stems
## drop secondary stem data
## drop records without a DBH
## drop records without spatial coordinates

goodData = dat$Status == 'alive' &
           dat$Stem == 'main' &
           !is.na(dat$DBH) &
           !is.na(dat$gx) &
           !is.na(dat$gy)  

ferp = dat[goodData,]

write.csv(ferp, file='./data/filtered_data/ferp_2007_filtered.csv', row.names=F)