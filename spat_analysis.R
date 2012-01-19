##Author: Dan McGlinn
##Purpose: to compare MaxEnt generate community to null models
##Outputs: the result of a null permutation routine.

library(danspkg)
library(vegan)
library(snowfall)

setwd('/home/danmcglinn/maxent')

source('spat_sim_vario_func.R')

clArgs = commandArgs(trailingOnly=TRUE)
if( length(clArgs) > 1){
  S = clArgs[1]
  N = clArgs[2]
  ncomm = clArgs[3]
  bisec = clArgs[4]
  transect = clArgs[5]
  dataType = clArgs[6]
  metricsToCalc = clArgs[7]
  direction = clArgs[8]
  tolerance = clArgs[9]
  name = clArgs[10]
}
if( !exists(as.character(substitute(S))) ){
  S = 10
  N = 1076
  ncomm = 200
  bisec = 13
  transect = FALSE
  dataType = 'both'
  metricsToCalc = 'all'
  direction = NA
  tolerance = NA
  name = NA
}

direction = ifelse(is.na(direction),'omnidirectional',as.numeric(direction))
tolerance = ifelse(is.na(tolerance),tolerance,as.numeric(tolerance)) 

if(is.na(name)){
  fileSuffix = ifelse(transect,
               paste('_S',S,'_N',N,'_C',ncomm,'_B',bisec,'_transect',sep=''),
               paste('_S',S,'_N',N,'_C',ncomm,'_B',bisec,'_grid',sep=''))
}
if(!is.na(name)){
  fileSuffix = ifelse(transect,
               paste('_',name,'_C',ncomm,'_B',bisec,'_transect',sep=''),
               paste('_',name,'_C',ncomm,'_B',bisec,'_grid',sep=''))
}

fileName = paste('simulated_comms',fileSuffix,'.txt',sep='')
comms = read.csv(file.path('./comms',fileName),header=T)
gc()

nperm = NULL
npar = NULL
RPargs = NULL
#nperm = 500
#npar = 8
#RPargs = c(TRUE,8,1e2,1e4,0.01,1)

dataType = ifelse(dataType=='both',c('binary','abu'),dataType)

sapply(dataType,function(x){
       calcMetrics(comms=comms,metricsToCalc=metricsToCalc,dataType=x,
                   direction=direction,tolerance=tolerance,nperm=nperm,npar=npar,
                   RPargs=RPargs,writeToFile=TRUE,fileSuffix=fileSuffix)})




