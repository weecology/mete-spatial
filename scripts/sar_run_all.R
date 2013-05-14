
setwd('~/maxent/spat')
dir.create('./data/')
dir.create('./sar')
dir.create('./figs')
dir.create('./scripts/log_files')

## ToDo: make it so that the sitenames are supplied once and 
## not within each individual script

setwd('./scripts')

## download and prepare data ---------------------------------------------
system('Rscript spat_download_empir_data.R', wait=TRUE)

## filter the empirical data 
system('Rscript spat_empir_data_filtering.R', wait=TRUE)

## calculate species-abundance distributions
system('Rscript spat_empir_sad.R', wait=TRUE)

## summarize the empirical data
system('Rscript spat_empir_data_summary.R', wait=TRUE)

## data analysis ---------------------------------------------------------

## analyze empirical observed SAR (approx 3 min)
system.time(system('Rscript spat_empir_sar.R', wait=TRUE))

## generate METE expected SAR (several hours to several days)
system('python spat_mete_sar.py >./log_files/mete_sar.log 2>&1', wait=FALSE)

## generate binomial expected SAR (approx 1.5 min)
system.time(system('Rscript spat_empir_expected_sars.R', wait=TRUE))

## aggregated results and generate figures ------------------------------
system('Rscript spat_sar_load_and_avg_data.R', wait=TRUE)

system('Rscript sar_figures.R', wait=TRUE)

