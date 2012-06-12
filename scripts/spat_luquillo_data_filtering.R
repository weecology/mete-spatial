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
