## Author: Dan McGlinn
## Purpose: to read in the .txt files of the BCI censues and output a cleaned
## .Rdata file that can then be used to perform calculations on. The filtering 
## performed in this script is for a biodiversity analysis and may not be
## appropriate for analyses targeted at other topics
## Metadata: bci50ha.doc 

setwd('~/datasets/CTFSplots/BCI')

fileNames = dir()

bciData = fileNames[grep('bci',fileNames)]
bciDataTxt = bciData[grep('txt',bciData)]


for(i in seq_along(bciDataTxt)){
  dat = read.table(bciDataTxt[i],sep='\t',header=TRUE)
  goodData = dat$Status=='alive' & !is.na(dat$DBH) & 
             !is.na(dat$gx) & !is.na(dat$gy) & 
             dat$Latin != 'Unidentified species' &
             dat$Stem != 'secondary'
  dat = dat[goodData ,]
  save(dat,file=paste(strsplit(bciDataTxt[i],'.txt'),'.Rdata',sep=''))
}

