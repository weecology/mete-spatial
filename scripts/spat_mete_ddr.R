## Purpose: to generate the analytical METE DDR predictions
## given the parameters in the empirically observed communities

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 0) {
  server = clArgs[1]
  site_index = clArgs[2]
}

shrtnames = read.table('./data/shrtnames.txt', colClasses='character')
grain_fine = as.numeric(read.table('./data/grain_fine.txt'))
bisect_fine = as.numeric(read.table('./data/bisect_fine.txt'))

i = site_index

Amin = 1
Amax = 2^bisect_fine[i]
shape = ifelse(log2(Amax) %% 2 == 0, 'sqr', 'rect')
abu_file = paste('./data/', shrtnames[i], '_sad.csv', sep='')
out_file = paste('./sorensen/', shrtnames[i], '_mete_sorensen.csv', sep='')
log_file = paste('./scripts/log_files/', shrtnames[i], '_mete_sor.log', sep='')
unit_distance = sqrt(grain_fine[i]) * n_pixels_wide(bisect_fine[i])
if (server == 'unc')
  cmd = paste('bsub -q day -o', log_file,
              'python ./scripts/spat_heap_ddr.py',
              Amin, Amax, shape, abu_file, out_file, 
              unit_distance, sep=' ')
else
  cmd = paste('python ./scripts/spat_heap_ddr.py',
              Amin, Amax, shape, abu_file, out_file, 
              unit_distance, '>', log_file, '2>&1', sep=' ')

system(cmd, wait=FALSE)


