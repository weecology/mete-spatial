## Author: Dan McGlinn
## Purpose: to read in the .txt files of the Sherman censues and output three 
## cleaned .csv files that can then be used to perform calculations on. The
## filtering performed in this script is for a biodiversity analysis and may not
## be appropriate for analyses targeted at other topics
## Metadata: sherman_files.doc
## There is a visual diagram of the Sherman plot in the following publication,
## Condit, R. et al. 2004. Tropical forest dynamics across a rainfall gradient
## and the impact of an El Nino dry season. Journal of Tropical Ecology, 20: 51-72.

setwd('/home/danmcglinn/CTFSplots/sherman')

dat = read.table('sherman.txt',header=TRUE)
spNames = read.table('shermansp.txt',header=TRUE)

goodData1 = dat$dbh1 > 0 & dat$recr1 == 'A' & dat$spcode != '*'
goodData2 = dat$dbh2 > 0 & dat$recr2 == 'A' & dat$spcode != '*'
goodData3 = dat$dbh3 > 0 & dat$recr3 == 'A' & dat$spcode != '*'

sherman1 = dat[goodData1,c('tag','spcode','x','y',
                          'dbh1','recr1','pom1','code1','mult1','date1')]
sherman2 = dat[goodData2,c('tag','spcode','x','y',
                          'dbh2','recr2','pom2','code2','mult2','date2')]
sherman3 = dat[goodData3,c('tag','spcode','x','y',
                          'dbh3','recr3','pom3','code3','mult3','date3')]
names(sherman1) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')
names(sherman2) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')
names(sherman3) = c('tag','spcode','x','y','dbh','recr','pom','code','mult','date')

write.csv(sherman1,file='sherman_census1_filtered.csv',row.names=FALSE)
write.csv(sherman2,file='sherman_census2_filtered.csv',row.names=FALSE)
write.csv(sherman3,file='sherman_census3_filtered.csv',row.names=FALSE)

