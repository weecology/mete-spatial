## Author: Dan McGlinn
## Purpose: to read in the .txt files of the BCI censues and output a cleaned
## csv file that can then be used to perform calculations on. The filtering 
## performed in this script is for a biodiversity analysis and may not be
## appropriate for analyses targeted at other topics
## Metadata: ~/datasets/CTFSplots/BCI/bci50ha.doc 

setwd('~/maxent/spat')

dat = read.table('./data/raw_data/bci_census7.txt', sep='\t', header=TRUE)

goodData = dat$Status=='alive' & !is.na(dat$DBH) & 
           !is.na(dat$gx) & !is.na(dat$gy) & 
           dat$Latin != 'Unidentified species' &
           dat$Stem != 'secondary'
dat = dat[goodData ,]

write.csv(dat, file='./data/filtered_data/bci_census7_filtered.csv', row.names=F)


