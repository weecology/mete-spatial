
setwd('~/maxent/spat')

dir.create('./comms')
dir.create('./sorensen')
dir.create('./sar')
dir.create('./scripts/log_files')

server = 'unc'

## empirical data analysis ------------------------------------------------

## filter the empirical data
system('Rscript spat_empir_data_filtering.R')

## summarize the empirical data
system('Rscript spat_empir_data_summary.R')

## simulate the empirical data
indices = paste(1:19, collapse=' ')
system(paste('Rscript spat_empir_comm_sim.R', server,
             "'", indices, "'", sep=" "))

## analyze emprical observed DDR 
npar = 8
nperm = 500
commName = as.character(read.table('../data/shrtnames.txt', colClasses='character'))
commName = paste("'", paste(commName, collapse=' '), "'", sep='')
system(paste("Rscript spat_empir_observed_ddr.R", server,
             npar, nperm, commName, sep=" "))

## analyze emprical simulated DDR
indices = paste("'", paste(1:19, collapse=' '), "'", sep='')
system(paste("Rscript spat_empir_simulated_ddr.R", server,
             indices, sep=" "))

## analyze empirical observed SAR 
system('Rscript spat_empir_sar.R')

## analyze empirical simulated SAR
sitename = as.character(read.table('../data/shrtnames.txt', colClasses='character'))
sitename = paste("'", sitename, "'", sep='')
system(paste('Rscript spat_mete_sar_avgs.R', sitename, sep=' '))

## generate binomial expected SAR 
system('Rscript spat_empir_expected_sars.R')

## generate METE expected SAR
system('chmod u+x spat_mete_sar.sh')
indices = paste("'", paste(1:19, collapse=' '), "'", sep='')
system(paste('./spat_mete_sar.sh', server, indices, sep=' '), wait=FALSE)

## parameter space analysis -----------------------------------------------
## simulate communities across parameter space
system('python gencomms.py')

## analyze parameter space communities DDR

## analyze parameter space communities SAR
system('Rscript spat_param_sar_avgs.R')
