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

write.csv(luqu, file='./data/filtered_data/luquillo_census4_filtered.csv',
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
## Purpose: to clean up the NC plots for analysis 
## Selecting the year prior to 1996 when Hurricane Fran hit
## the recent disturbances were Umstead tornado 1988, Hurricane Hugo 1989, and 
## Hurricane Fran 1996
## Metadata: is in the raw data files and in ~/datasets/NC_LTREB/Maps/Readme_maps.txt 
## C = CC** = Condition code for hear of preceding diameter D**
## 1=live, OK; 2=dead; 3=missing; 4=died back below breast ht; 5=cut
## 6=damaged by Hurrricane Fran in 97; consult Hurricane codes
## coordinates are sometimes in decimeters and somtimes in meters, always check
## the readme file.

setwd('~/maxent/spat')

dat_names = paste('m', c('04','07','12','13',91:94), sep='')
plc_names = c('graveyard','landsend','rocky','bormann',
              'woodbridge','baldmnt','bryan','bigoak')

file_names = paste(dat_names, plc_names, sep='_')
skip_lines = c(79, 20, 58, 49, 10, 10, 9, 87)

dat = vector('list', length(dat_names))
for (i in seq_along(dat_names)) {  
  file_path = paste('./data/raw_data/', file_names[i], '.csv', sep='')
  tmp = read.csv(file_path, skip=skip_lines[i] + 1, header=F, na.strings='.')
  names(tmp) = as.matrix(read.csv(file_path, skip=skip_lines[i], nrows=1, header=F))
  dat[[i]] = tmp
}
names(dat) = dat_names

head(dat[[1]])  ## c(78, 82, 89, 92, 97)
head(dat[[2]])  ## c(78, 82, 89, 93, 97)
head(dat[[3]])  ## c(78, 82, 85, 90, 97)
head(dat[[4]])  ## c(50, 74, 82, 89, 93, 97, 00)
head(dat[[5]])  ## c(84, 86, 91, 98)
head(dat[[6]])  ## c(84, 86, 91)
head(dat[[7]])  ## c(86, 91, 98)
head(dat[[8]])  ## c(86, 90, 93, 97, 01)

yrs = c(92, 93, 90, 93, 91, 91, 91, 93)

dat_filter = dat

i = 1
tmp = dat_filter[[i]]

head(tmp)
true = !is.na(tmp$D92) & 
       tmp[,13] == 1 &
       tmp$X <= 1000


tmp$X = tmp$X / 10
tmp$Y = tmp$Y / 10
plot(tmp$X[true], tmp$Y[true])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D92')]
###
i = 2

tmp = dat_filter[[i]]

head(tmp)
true = !is.na(tmp$D93) & 
       tmp[,12] == 1 &
       tmp$X <= 1300 & 
       tmp$Y <= 650 

tmp$X = tmp$X / 10
tmp$Y = tmp$Y / 10
plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D93')]
###
i = 3

tmp = dat_filter[[i]]

head(tmp)
true = !is.na(tmp$D90) & 
       tmp[,13] == 1 &
       tmp$X <= 1200 & 
       tmp$Y <= 1200

tmp$X = tmp$X / 10
tmp$Y = tmp$Y / 10
plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D90')]
###
i = 4

tmp = dat_filter[[i]]

head(tmp)
names(tmp)
true = !is.na(tmp$D93) &
       tmp[,13] == 1

plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D93')]
###
i = 5

tmp = dat_filter[[i]]

head(tmp)
names(tmp)
true = !is.na(tmp$D91) & 
       tmp[,12] == 1 & 
       tmp$X <= 710 &
       tmp$Y <= 710

tmp$X = tmp$X / 10
tmp$Y = tmp$Y / 10
plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D91')]
###
i = 6

tmp = dat_filter[[i]]

head(tmp)
names(tmp)
true = !is.na(tmp$D91) & 
       tmp[,12] == 1 
  
tmp$X = tmp$X / 10
tmp$Y = tmp$Y / 10
plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D91')]
###
i = 7

tmp = dat_filter[[i]]

head(tmp)
names(tmp)
true = !is.na(tmp$D91) & 
       tmp$CC91 == 1 & 
       tmp$X <= 185 &
       tmp$Y <= 185 / 2

plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D91')]
###
i = 8

tmp = dat_filter[[i]]

head(tmp)
names(tmp)
true = !is.na(tmp$D93) & 
       tmp[,10] == 1 & 
       tmp$X <= 200 

plot(tmp$X[true], tmp$Y[true])

head(tmp[true,])

dat_filter[[i]] = tmp[true, c('ID','SPEC','X','Y','D93')]

###
## clean up taxonomy
sp_lookup = read.csv('./data/raw_data/species.csv')

## place all datasets into a single flat file

tmp = dat_filter
names(tmp[[1]]) = c("ID", "SPEC", "X", "Y", "DBH")
for (i in seq_along(tmp)) {
  if (i == 1)
    dat_flat = data.frame(dat=dat_names[i], loc = plc_names[i], tmp[[i]])
  else {
    names(tmp[[i]]) = names(tmp[[1]])
    dat_flat = rbind(dat_flat, data.frame(dat=dat_names[i], loc = plc_names[i],
                                          tmp[[i]]))
  }  
}

index = match(dat_flat$SPEC, sp_lookup$SPEC)
SciName = ifelse(is.na(index), NA, as.character(sp_lookup$ScientificName[index]))
dat_flat = data.frame(dat_flat, ScientificName = SciName)

head(dat_flat)

maps_name_errors = dat_flat[is.na(dat_flat$Sci), ]
write.csv(maps_name_errors, file='./data/maps_name_errors.csv', row.names=F)

## bring in information from Bob Peet and Michael Lee 

correct_key = rbind(c("0OFL", "COFL"), 
                    c("CHIV", "CHVI"),
                    c("CPF;", "COFL"),
                    c("EUON", "EUAM"),  
                    c("ILAX", "ILEX"),
                    c("JPRS", "PRSE"),
                    c("MACO", "VIPR"),
                    c("PLEC", "PIEC"))

colnames(correct_key) = c('SPEC','SPEC_corrected')

dat_flat$SPEC = as.character(dat_flat$SPEC)
dat_flat = dat_flat[ , -8]
for (i in 1:nrow(correct_key)) 
  dat_flat$SPEC[dat_flat$SPEC == correct_key[i, 1]] = correct_key[i, 2]

index = match(dat_flat$SPEC, sp_lookup$SPEC)
SciName = ifelse(is.na(index), NA, as.character(sp_lookup$ScientificName[index]))
dat_flat = data.frame(dat_flat, ScientificName = SciName)

maps_name_errors = dat_flat[is.na(dat_flat$Sci), ]
## drop all stems in the maps_name_errors file
false = dat_flat$SPEC %in% unique(maps_name_errors$SPEC)
dat_flat = dat_flat[!false, ]

## now sort out what to do with the records identified only to genus

sp_dat = dat_flat[grep(" sp.", dat_flat$Sci),]
genera = sub(' sp.','', unique(sp_dat$Sci)) 

ge = genera[1]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## so Fraxinus sp. is the only identifier for Ash trees in the entire dataset

ge = genera[2]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## in all the datsets except m13 it looks like just dropping the Carya sp. 
## records is the best approach

ge = genera[3]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## the Pinus sp. can be dropped

ge = genera[4]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## the Quercus sp. can be dropped

ge = genera[5]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## Crataegus sp. is the only identifier for Hawthorns in the dataset

ge = genera[6]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## Vaccinium sp. is the only blueberry plant id'ed in that particular dataset

ge = genera[7]
ge_dat = dat_flat[grep(ge, dat_flat$Sci),]
x = tapply(ge_dat$Sci, list(ge_dat$dat), table)
lapply(x, function(y) y[y>0])
## Ilex sp. can be dropped

## So the conclusions are:
## 1) keep Fraxinus sp.
## 2) drop Carya sp. in all datasets except m13 where all Caryas should be
##    lumped to Carya sp.
ge = 'Carya'
index = grep(ge, dat_flat$Sci)
index = index[dat_flat$dat[index] == "m13"]
dat_flat$SPEC[index] = "CARY"
##
ge = 'Carya sp.'
index = grep(ge, dat_flat$Sci)
index = index[dat_flat$dat[index] != "m13"]
dat_flat = dat_flat[-index, ]
## 3) drop Pinus sp.
dat_flat = dat_flat[dat_flat$SPEC != "PINU", ]
## 4) drop Quercus sp.
dat_flat = dat_flat[dat_flat$SPEC != "QUER", ]
## 5) keep Crategus sp.
## 6) keep Vaccinium sp.
## 7) drop Ilex sp. 
dat_flat = dat_flat[dat_flat$SPEC != "ILEX", ]

## check on the clean up
sp_dat = dat_flat[grep(" sp.", dat_flat$Sci),]
table(sp_dat$SPEC)
## this looks correct

## now output dat_flat into seperate files for each dataset

## output filtered files
for (i in seq_along(dat_names)) {
  tmp = dat_flat[dat_flat$dat == dat_names[i], -c(1, 2, 8)]
  write.csv(tmp, file= paste('./data/filtered_data/', file_names[i],'_19', 
                             yrs[i], '_filtered.csv', sep=''), row.names=FALSE)
}


