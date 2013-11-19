#!/usr/bin/Rscript

## install custom python dependencies

system('git clone git@github.com:weecology/METE.git')
setwd('./METE')
system('python setup.py install')
setwd('..')

system('git clone git@github.com:weecology/macroecotools.git')
setwd('./macroecotools')
system('python setup.py install')
setwd('..')

## Additional python packages required include:
## matplotlib, mpmath, numpy, scipy

## install R dependencies
install.packages(c('vegan', 'RCurl', 'bigmemory'))

dir.create('./comms')
dir.create('./sorensen')
dir.create('./scripts/log_files')

server = 'unc'

## download and prepare data ---------------------------------------------
system('Rscript ./scripts/spat_download_empir_data.R')

## get names of datasets
system('Rscript ./scripts/spat_get_sitenames.R')

## filter the empirical data 
system('Rscript ./scripts/spat_empir_data_filtering.R')

## create multiscale site x species matrics (several minutes to several hours)
system('Rscript ./scripts/spat_empir_gen_comm_mat.R')

## calculate empirical SAD
system('Rscript ./scripts/spat_empir_sad.R')

## summarize the empirical data
system('Rscript ./scripts/spat_empir_data_summary.R')

## simulate communities that are parameterized by the empirical communities (slow)
indices = paste(1:19, collapse=' ')
system(paste('Rscript spat_empir_comm_sim.R', server,
             "'", indices, "'", sep=" "), wait=FALSE)

## data analysis ---------------------------------------------------------

## analyze emprical observed DDR (very slow)
## 1) analyze the multivariate DDR with the permutation based RP for abundance data
npar = 8
nperm = 500
commName = 'all'
method = 'multi'
dataType = 'both'
system(paste("Rscript spat_empir_observed_ddr.R", server,
             npar, nperm, commName, method, dataType, sep=" "),
       wait=FALSE)

## 2) analyze the univariate (i.e, species specific) DDR just for the abundance data.  
## Note: the number of processers the analysis should run over will depend on
## availability of system memory and the runs may become very bogged down if
## care is not taken
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

## compute the empirical observed DDR using spat_empir_analysis.R
system('Rscript ./scripts/spat_empir_observed_ddr.R')

## aggregate, save, and graph the simulated and empirical results
system('Rscript ./scripts/spat_sim_results_summary.R')
system('Rscript ./scripts/spat_empir_results_summary.R')

## compare the fit of the simulated ddr to the empirical ddr
system('Rscript ./scripts/spat_ddr_comparision.R')

## generate figures ------------------------------------------------------

system('Rscript ./scripts/ddr_figures.R')
