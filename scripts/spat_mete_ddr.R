## Purpose: to generate the analytical METE DDR predictions
## given the parameters in the empirically observed communities

print('Submitting jobs to compute the analytical METE DDRs, ...')

source('./scripts/spat_functions.R')

clArgs = commandArgs(trailingOnly=TRUE)
if (length(clArgs) > 0) {
  server = clArgs[1]
  site_index = clArgs[2]
  sadType = clArgs[3]
}

if (sadType == 'both')
  sadType = c('empirSAD','meteSAD')

shrtnames = read.table('./data/shrtnames.txt', colClasses='character')
grain_fine = as.numeric(read.table('./data/grain_fine.txt'))
bisect_fine = as.numeric(read.table('./data/bisect_fine.txt'))
bisect_coarse = as.numeric(read.table('./data/bisect_coarse.txt'))

if (site_index != 'all')
  site_index = as.integer(site_index)

if (site_index == 'all') { 
  site_index = 1:length(shrtnames)
}

nice = FALSE # make TRUE to use nice when submitting jobs
wait = TRUE  # make FALSE to submit all the jobs at once

for(i in site_index) {
  for(j in sadType) {
    unit_distance = sqrt(grain_fine[i]) * n_pixels_wide(bisect_fine[i])
    abu_file = paste('./data/', shrtnames[i], '_sad.csv', sep='')    
    if (j == 'meteSAD') {
      out_file = paste('./sorensen/', shrtnames[i], 
                       '_mete_sor.csv', sep='')
      log_file = paste('./scripts/log_files/', shrtnames[i],
                       '_mete_sor.log', sep='')
    }
    if (j == 'empirSAD') {
      out_file = paste('./sorensen/', shrtnames[i], 
                       '_empirSAD_mete_sor.csv', sep='')
      log_file = paste('./scripts/log_files/', shrtnames[i],
                       '_empirSAD_mete_sor.log', sep='')    
    }    
    if (server == 'LSF') {
      cmd = paste('bsub -q week -M 10 -o', log_file, '-J', shrtnames[i],
                  'python ./scripts/spat_heap_ddr.py',
                  bisect_fine[i], bisect_coarse[i], 
                  j, abu_file, out_file, 
                  unit_distance, sep=' ')
    }
    else {
      if (nice)
        cmd = 'nice python'
      else
        cmd = 'python'
      cmd = paste(cmd, './scripts/spat_heap_ddr.py',
                  bisect_fine[i], bisect_coarse[i],
                  j, abu_file, out_file, 
                  unit_distance, '>', log_file, '2>&1', sep=' ')
    
    system(cmd, wait=wait)
  }
}

print('Submitting jobs to compute the analytical METE DDRs, complete!')
