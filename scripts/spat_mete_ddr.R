## Purpose: to generate the analytical METE DDR predictions
## given the parameters in the empirically observed communities

setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 0) {
  server = clArgs[1]
  site_index = as.integer(clArgs[2])
  sadType = clArgs[3]
}

shrtnames = read.table('./data/shrtnames.txt', colClasses='character')
grain_fine = as.numeric(read.table('./data/grain_fine.txt'))
bisect_fine = as.numeric(read.table('./data/bisect_fine.txt'))
bisect_coarse = as.numeric(read.table('./data/bisect_coarse.txt'))

i = site_index

unit_distance = sqrt(grain_fine[i]) * n_pixels_wide(bisect_fine[i])
abu_file = paste('./data/', shrtnames[i], '_sad.csv', sep='')
log_file = paste('./scripts/log_files/', shrtnames[i], '_mete_sor.log', sep='')

if (sadType == 'meteSAD')
  out_file = paste('./sorensen/', shrtnames[i], 
                   '_mete_sor.csv', sep='')
if (sadType == 'empirSAD')
  out_file = paste('./sorensen/', shrtnames[i], 
                   '_empirSAD_mete_sor.csv', sep='')

if (server == 'unc')
  cmd = paste('bsub -q week -x -o', log_file, '-J', shrtnames[i],
              'python ./scripts/spat_heap_ddr.py',
              bisect_fine[i], bisect_coarse[i], 
              sadType, abu_file, out_file, 
              unit_distance, sep=' ')
if (server == 'usu')
  cmd = paste('nice python ./scripts/spat_heap_ddr.py',
              bisect_fine[i], bisect_coarse[i],
              sadType, abu_file, out_file, 
              unit_distance, '>', log_file, '2>&1', sep=' ')

system(cmd, wait=FALSE)


