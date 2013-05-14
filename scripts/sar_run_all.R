
setwd('~/maxent/spat')
dir.create('./data/')
dir.create('./sar')
dir.create('./figs')

## ToDo: make it so that the sitenames are supplied once and 
## not within each individual script

setwd('./scripts')

## download and prepare data ---------------------------------------------
## requires R package: RCurl
system('Rscript spat_download_empir_data.R')

## filter the empirical data 
system('Rscript spat_empir_data_filtering.R')

## calculate empirical SAD
system('Rscript spat_empir_sad.R')

## summarize the empirical data
system('Rscript spat_empir_data_summary.R')

## data analysis ---------------------------------------------------------

## compare METE and empirical SAD
system('python spat_obs_pred_sad.py')

## analyze empirical observed SAR (approx 3 min)
system('Rscript spat_empir_sar.R')

## generate METE expected SAR (several hours to several days)
system('python spat_mete_sar.py')

## generate binomial expected SAR (approx 1.5 min)
system('Rscript spat_empir_expected_sars.R')

## aggregated results and generate figures ------------------------------
system('Rscript spat_sar_load_and_avg_data.R')

system('Rscript sar_figures.R')

system('python spat_plot_obs_pred_sad.py')
