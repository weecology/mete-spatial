## Author: Dan McGlinn
## Purpose: to read in the .txt files of the cocoli censues and output three 
## cleaned .csv files that can then be used to perform calculations on. The
## filtering performed in this script is for a biodiversity analysis and may not
## be appropriate for analyses targeted at other topics
## Metadata: cocoli_files.doc
## There is a visual diagram of the Sherman plot in the following publication,
## which is fairly close to the cocoli plot except for two 40 m strips 
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Niño dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('/home/danmcglinn/CTFSplots/cocoli')

dat = read.table('cocoli.txt',header=TRUE)
spNames = read.table('cocolisp.txt',header=TRUE)

goodData1 = dat$dbh1 > 0 & dat$recr1 == 'A' & dat$spcode != '*'
goodData2 = dat$dbh2 > 0 & dat$recr2 == 'A' & dat$spcode != '*'
goodData3 = dat$dbh3 > 0 & dat$recr3 == 'A' & dat$spcode != '*'

cocoli1 = dat[goodData1,c('tag','spcode','x','y',
                          'dbh1','recr1','pom1','code1','mult1','date1')]
cocoli2 = dat[goodData2,c('tag','spcode','x','y',
                          'dbh2','recr2','pom2','code2','mult2','date2')]
cocoli3 = dat[goodData3,c('tag','spcode','x','y',
                          'dbh3','recr3','pom3','code3','mult3','date3')]
names(cocoli1) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')
names(cocoli2) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')
names(cocoli3) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')

write.csv(cocoli1,file='cocoli_census1_filtered.csv',row.names=FALSE)
write.csv(cocoli2,file='cocoli_census2_filtered.csv',row.names=FALSE)
write.csv(cocoli3,file='cocoli_census3_filtered.csv',row.names=FALSE)

