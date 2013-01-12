## Purpose: to compute the empirical simulated DDR 

setwd('~/maxent/spat/scripts')

clArgs = commandArgs(trailingOnly=TRUE)

server = clArgs[1]
indices = clArgs[2]

## arguments for job
S = read.table('../data/S_vals.txt')
N = read.table('../data/N_vals.txt')
bisect_fine = read.table('../data/bisect_fine.txt')
bisect_coarse = read.table('../data/bisect_coarse.txt')
grain_fine = read.table('../data/grain_fine.txt')
names = as.character(read.table('../data/shrtnames.txt', colClasses='character'))

sadType = c("meteSAD", "empirSAD")
dataType = c("abu", "binary")
metrics = "sorensen"

for (i in indices) {
  for (j in sadType) {
    for (k in dataType) {
      for (m in metrics) {
        log_file = paste('./log_files/error_sim_analysis_', names[i],
                         '_', j, '.log', sep='')
        if (j == 'meteSAD')
          name = names[i]
        if (j == 'empirSAD')
          name = paste(names[i], '_empirSAD', sep='')
        if ( sever == unc) 
          system(paste('bsub -q week -M 8 -J', names[i],
                       '-o', log_file, 'Rscript spat_analysis.R',
                       S[i], N[i], 200, bisect_fine[i], bisect_coarse[i],
                       grain_fine[i], FALSE, k, m, NA, NA,
                       name, TRUE, sep=' '))
        else
          system(paste('Rscript spat_analysis.R', S[i], N[i], 200,
                       bisect_fine[i], bisect_coarse[i],
                       grain_fine[i], FALSE, k, m, NA, NA,
                       name, TRUE, '>', log_file, '2>&1 &', sep=' '))
      }    
    }
  }
}
