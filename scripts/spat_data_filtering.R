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

##------------------------------------------------------------------------------

## Purpose: to read in the .txt files of the cocoli censues and output a 
## cleaned .csv files that can then be used to perform calculations on. The
## filtering performed in this script is for a biodiversity analysis and may not
## be appropriate for analyses targeted at other topics
## Metadata: cocoli_files.doc
## There is a visual diagram of the Sherman plot in the following publication,
## which is fairly close to the Cocoli plot except for two 40 m strips 
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('~/maxent/spat')

dat = read.table('./data/raw_data/cocoli.txt', header=TRUE)
spNames = read.table('.data/raw_data/cocolisp.txt', header=TRUE)

goodData = dat$dbh3 > 0 & dat$recr3 == 'A' & dat$spcode != '*'

cocoli = dat[goodData, c('tag','spcode','x','y',
                        'dbh3','recr3','pom3','code3','mult3','date3')]
names(cocoli) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')

write.csv(cocoli, file='./data/filtered_data/cocoli_census3_filtered.csv',
          row.names=FALSE)

##------------------------------------------------------------------------------

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

##------------------------------------------------------------------------------

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

##------------------------------------------------------------------------------

## Purpose: to clean up the 2007 Ferp woodly plant dataset for analysis of 
## spatial biodiversity patterns. 
## Metadata: http://ferp.ucsc.edu/ferp/about/using/

setwd('~/maxent/spat')

dat= read.delim('./data/raw_data/FERP_CTFS_2007_data.txt')

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

range(ferp$gx)
range(ferp$gy)
## drop records that are less than zero
ferp = ferp[ferp$gx > 0 & ferp$gy > 0, ]

write.csv(ferp, file='./data/filtered_data/ferp_2007_filtered.csv', row.names=F)

##------------------------------------------------------------------------------

## Purpose: to clearn the Luquillo forest dataset for analysis of 
## spatial biodiversity patterns, this survey is from 2006/2007
## Metadata: http://luq.lternet.edu/data/luqmetadata119
## see the following page for the status codes defined:
## http://luq.lternet.edu/data/variable/codes

setwd('~/maxent/spat')

dat1 = read.csv('./data/raw_data/LFDP_Census4-Part1.txt')
dat2 = read.csv('./data/raw_data/LFDP_Census4-Part2.txt')
dat = rbind(dat1, dat2)

## drop secondary stem data
## drop unidentified species
## drop records without spatial coordinates
## drop records without a DBH
goodData = dat$TAG == dat$STEMTAG &
           dat$SPCODE != "" &
           dat$SPCODE != "???" &
           !is.na(dat$GX) &
           !is.na(dat$GY) & 
           !is.na(dat$DBH)

luqu = dat[goodData, ]

range(luqu$GX)
# the max spatial x should be 320 not 320.02
# change this record 
luqu$GX[luqu$GX == max(luqu$GX)] = 319.99

write.csv(luqu, file='./data/filtered_data/luqu_census4_filtered.csv',
          row.names=F)

##------------------------------------------------------------------------------

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

##------------------------------------------------------------------------------


