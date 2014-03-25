#!/usr/bin/Rscript

print('Analysis Begin')
start.time = proc.time()

## Additional python packages required include:
## matplotlib, mpmath, numpy, scipy

## Additional R packages required include:
# install.packages(c('vegan', 'RCurl', 'bigmemory', 'snowfall', 'rlecuyer'))

dir.create('./comms')
dir.create('./sorensen')
dir.create('./scripts/log_files')
dir.create('./figs')

## specify server type to run code on
## options are 'local' or 'LSF'
## local = local machine 
## LSF = load sharing facility
server = 'local' 
                 
## specify analysis type
## options are 'quick' or 'full'
## quick --> provides example results for the oosting and ucsc sites
## full --> provies all results for all 16 sites Warning: takes a long time and requires data not publicly provided
analysis_type = 'quick' 

set_analysis_params = function(analysis_type) {
  params = list()
  if (analysis_type == 'quick') {
    params$sites = c('oosting', 'ucsc')
    params$ncomm = 10
    params$nperm = 20
    params$npar = 2 
  }
  if (analysis_type == 'full') {
    params$sites = NULL
    params$ncomm = 200
    params$nperm = 500
    params$npar = 4
  }
  return(params)
}

params = set_analysis_params(analysis_type)

## download and prepare data ---------------------------------------------
system(paste('Rscript ./scripts/spat_download_empir_data.R', 
             paste(params$sites, collapse=' ')))

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
sites = read.table('./data/shrtnames.txt', colClasses='character')[1, ]
sites = paste(as.character(sites), collapse=' ')
#number of community matrices to simulatie
ncomm = params$ncomm
system(paste('Rscript ./scripts/spat_empir_comm_sim.R', server, ncomm, sites),
       wait=TRUE)

## data analysis ---------------------------------------------------------

## analyze the empirical multivariate DDR with the individual-based
## RPM for abundance data
## Warning! This calculation requires tremendous computataion resources.
## Take care before submitting job that enough ram and processors are available.
## If run in serial this job will likely take 2 weeks to complete.
npar = params$npar   # number of processors to run code on.
nperm = params$nperm # specify as NA if RPM is not to be generated
commName = 'all'
method = 'multi'
dataType = 'both'
system(paste('Rscript ./scripts/spat_empir_observed_ddr.R', server,
             npar, nperm, commName, method, dataType, sep=' '),
       wait=TRUE)

## analyze the METE simulated DDRs (also a time intensive calculation)
## Warning! this script submits many jobs to a server simleantously
## care must be taken here not to crash your machine if a job
## schedular is not being implemented
## Warning! this will take at a minimum several days to complete
commName = 'all'
dataType = 'both'
bisect = TRUE
ncomm = params$ncomm
system(paste("Rscript ./scripts/spat_empir_simulated_ddr.R", server,
             commName, ncomm, dataType, bisect, sep=" "), wait=TRUE)

## generate analytical mete DDR expectation for the empirical communities---------
## this is extremely time intensive and may fail for large study 
## sites due to the large amounts of RAM that are needed
gen_ddr_analytical = FALSE
if (gen_ddr_analytical) {
  site_index = 'all'
  sadType = 'both'
  system(paste('Rscript ./scripts/spat_mete_ddr.R', server, site_index, sadType),
         wait=TRUE)
}

## aggregate and save the simulated and empirical results
system('Rscript ./scripts/spat_empir_results_summary.R')
system('Rscript ./scripts/spat_sim_results_summary.R')

## compare the fit of the simulated ddr to the empirical ddr
system('Rscript ./scripts/spat_ddr_comparison.R')

## generate figures ------------------------------------------------------------
system('Rscript ./scripts/ddr_figures.R')

print('Figures are stored in the directory ./figs/')

end.time = proc.time()
print ('Analysis Complete')
print (end.time - start.time)

