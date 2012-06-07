## Author: Dan McGlinn
## Purpose: to read in the .csv file of the crosstimbers 1998 censues and output a
## cleaned .csv file that can then be used to perform calculations on. The
## filtering performed in this script is for a biodiversity analysis and may not
## be appropriate for analyses targeted at other topics
## Metadata: CrosstimbersMasterDataSheet1-20-12.xlsx

setwd('~/maxent/spat')

dat = read.csv('./data/raw_data/crosstimbers1998.csv')

head(dat)
names(dat) = c('T','Q','N','sp','x','y','dbh','comments')
sort(unique(dat$T))
dat$T[as.character(dat$T)=='I'] = 'l'
sort(unique(dat$sp))
dat$sp[as.character(dat$sp)=='QUST '] = 'QUST'
dat = dat[,-ncol(dat)]
write.csv(dat,file='./data/filtered_data/cross1998_filtered.csv',row.names=FALSE)


