
setwd('~/maxent/spat')

dir.create('./comms')
dir.create('./sorensen')
dir.create('./sar')
dir.create('./scripts/log_files')

server = 'unc'

## empirical data analysis ------------------------------------------------

## filter the empirical data (quick)
system('Rscript spat_empir_data_filtering.R', wait=TRUE)

## summarize the empirical data (quick)
system('Rscript spat_empir_data_summary.R', wait=TRUE)

## simulate the empirical data (slow)
indices = paste(1:19, collapse=' ')
system(paste('Rscript spat_empir_comm_sim.R', server,
             "'", indices, "'", sep=" "), wait=FALSE)

## analyze emprical observed DDR (very slow)
npar = 8
nperm = 500
commName = as.character(read.table('../data/shrtnames.txt', colClasses='character'))
commName = paste("'", paste(commName, collapse=' '), "'", sep='')
system(paste("Rscript spat_empir_observed_ddr.R", server,
             npar, nperm, commName, sep=" "), wait=FALSE)

## analyze emprical simulated DDR (very slow)
indices = paste("'", paste(1:19, collapse=' '), "'", sep='')
system(paste("Rscript spat_empir_simulated_ddr.R", server,
             indices, sep=" "), wait=FALSE)

## generate analytical mete DDR expectation for the empirical communities

system('Rscript spat_mete_ddr.R', wait=FALSE)

## analyze empirical observed SAR (3 min)
system('Rscript spat_empir_sar.R', wait=FALSE)

## analyze empirical simulated SAR 
sitename = as.character(read.table('../data/shrtnames.txt', colClasses='character'))
for (i in sitename)
  system(paste('Rscript spat_mete_sar_avgs.R', i, sep=' '),
         wait=FALSE)

## generate binomial expected SAR 
system('Rscript spat_empir_expected_sars.R', wait=FALSE)

## generate METE expected SAR
system('chmod u+x spat_mete_sar.sh')
indices = paste("'", paste(1:19, collapse=' '), "'", sep='')
system(paste('./spat_mete_sar.sh', server, indices, sep=' '), wait=FALSE)
## need to update spat_mete_sar.py so that it computes the results
## for both the inferred and empirical SAD


## parameter space analysis -----------------------------------------------
## simulate communities across parameter space
system('python gencomms.py', wait=FALSE)

## analyze parameter space communities DDR

## analyze parameter space communities SAR
system('Rscript spat_param_sar_avgs.R', wait=FALSE)
