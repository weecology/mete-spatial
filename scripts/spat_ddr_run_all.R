
setwd('~/maxent/spat')

dir.create('./comms')
dir.create('./sorensen')
dir.create('./scripts/log_files')

server = 'unc'

## empirical data analysis ------------------------------------------------
setwd('./scripts')

## filter the empirical data (quick)
system('Rscript spat_empir_data_filtering.R', wait=TRUE)

## summarize the empirical data (quick)
system('Rscript spat_empir_data_summary.R', wait=TRUE)

## simulate the empirical data (slow)
indices = paste(1:19, collapse=' ')
system(paste('Rscript spat_empir_comm_sim.R', server,
             "'", indices, "'", sep=" "), wait=FALSE)

## analyze emprical observed DDR (very slow)
## 1) analyze the multivariate DDR with 
## the permutation based RP for abundance data
npar = 8
nperm = 500
commName = 'all'
method = 'multi'
dataType = 'both'
system(paste("Rscript spat_empir_observed_ddr.R", server,
             npar, nperm, commName, method, dataType, sep=" "),
       wait=FALSE)

## 2) analyze the univariate (i.e, species specific) DDR 
## just for the abundance data.  
## Note: the number of processers the analysis should run over
## will depend on availability of system memory and the runs may
## become very bogged down if care is not taken
npar = 5
nperm = NA
commName = 'all'
method = 'uni'
dataType = 'abu'
system(paste("Rscript spat_empir_observed_ddr.R", server,
             npar, nperm, commName, method, dataType, sep=" "),
       wait=FALSE)

## generate analytical mete DDR expectation for the empirical communities

server = 'usu'
site_index = 'all'
sadType = 'both'
system(paste('Rscript spat_mete_ddr.R', server, site_index, sadType), wait=FALSE)

## analyze emprical simulated DDR (very slow)
commName = 'all'
dataType = 'both'
bisect = TRUE

system(paste("Rscript spat_empir_simulated_ddr.R", server,
             commName, dataType, bisect, sep=" "), wait=FALSE)

