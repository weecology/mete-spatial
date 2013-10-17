#!/usr/bin/Rscript

## install dependencies
install.packages('RCurl')

system('git clone git@github.com:weecology/METE.git')
setwd('./METE')
system('python setup.py install')
setwd('..')

## download and prepare data ---------------------------------------------
## requires R package: RCurl
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

## data analysis ---------------------------------------------------------

## compare METE and empirical SAD
system('python ./scripts/spat_obs_pred_sad.py')

## analyze empirical observed SAR (approx 3 min)
system('Rscript ./scripts/spat_empir_sar.R')

## generate METE expected SAR (several minutes to several days)
system('python ./scripts/spat_mete_sar.py')

## generate binomial expected SAR (approx 1.5 min)
system('Rscript ./scripts/spat_empir_expected_sars.R')

## aggregated results and generate figures ------------------------------
system('Rscript ./scripts/spat_sar_load_and_avg_data.R')

system('Rscript ./scripts/sar_figures.R')

#system('python ./scripts/spat_plot_obs_pred_sad.py')

