#!/usr/bin/Rscript

get_dependencies = function(pkgs, ...) {
  ## arguments
  ## pkgs : a character vector of the names of the packages to check for
  ## ... : other optiomal arguments to supply to install.packages()
  local_pkgs = installed.packages()[ , 1]
  for (p in pkgs) {
    if (is.element(p, local_pkgs)) 
      print(paste(p, 'already installed'))
    else
      install.packages(p, ...)
  }
} 

get_dependencies('RCurl')

setwd('~/maxent/spat/scripts')

## download and prepare data ---------------------------------------------
## requires R package: RCurl
system('Rscript spat_download_empir_data.R')

## get names of datasets
system('Rscript spat_get_sitenames.R')

## filter the empirical data 
system('Rscript spat_empir_data_filtering.R')

## create multiscale site x species matrics (several minutes to several hours)
system('Rscript spat_empir_gen_comm_mat.R')

## calculate empirical SAD
system('Rscript spat_empir_sad.R')

## summarize the empirical data
system('Rscript spat_empir_data_summary.R')

## data analysis ---------------------------------------------------------

## compare METE and empirical SAD
system('python spat_obs_pred_sad.py')

## analyze empirical observed SAR (approx 3 min)
system('Rscript spat_empir_sar.R')

## generate METE expected SAR (several minutes to several days)
system('python spat_mete_sar.py')

## generate binomial expected SAR (approx 1.5 min)
system('Rscript spat_empir_expected_sars.R')

## aggregated results and generate figures ------------------------------
system('Rscript spat_sar_load_and_avg_data.R')

system('Rscript sar_figures.R')

system('python spat_plot_obs_pred_sad.py')

