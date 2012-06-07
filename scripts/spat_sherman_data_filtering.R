## Author: Dan McGlinn
## Purpose: to read in the .txt files of the Sherman censues and output a 
## cleaned .csv file that can then be used to perform calculations on. The
## filtering performed in this script is for a biodiversity analysis and may not
## be appropriate for analyses targeted at other topics
## Metadata: sherman_files.doc
## There is a visual diagram of the Sherman plot in the following publication,
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('~/maxent/spat')

dat = read.table('./data/raw_data/sherman.txt', header=TRUE)
spNames = read.table('./data/raw_data/shermansp.txt', header=TRUE)

## only use the data from the 3rd census
goodData = dat$dbh3 > 0 & dat$recr3 == 'A' & dat$spcode != '*'

sherman = dat[goodData, c('tag','spcode','x','y',
                         'dbh3','recr3','pom3','code3','mult3','date3')]
names(sherman) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')

write.csv(sherman, file='./data/filtered_data/sherman_census3_filtered.csv',
          row.names=FALSE)

