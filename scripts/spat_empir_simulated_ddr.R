## Purpose: to compute the empirical simulated DDR 

setwd('~/maxent/spat/scripts')

clArgs = commandArgs(trailingOnly=TRUE)

server = clArgs[1]
commName = clArgs[2]
dataType = clArgs[3]
bisect = clArgs[4]
memory = clArgs[5]
iteration = clArgs[6]

if (dataType == 'both') {
  dataType = c("abu", "binary")
}
if (is.na(memory)) {
  memory = 10
}
if (is.na(iteration)) {
  iteration = 1
}

## arguments for job
S = read.table('../data/S_vals.txt')
N = read.table('../data/N_vals.txt')
bisect_fine = read.table('../data/bisect_fine.txt')
bisect_coarse = read.table('../data/bisect_coarse.txt')
grain_fine = read.table('../data/grain_fine.txt')
names = as.character(read.table('../data/shrtnames.txt', colClasses='character'))

commName = unlist(strsplit(commName, ' '))
if (commName[1] == 'all') {
  commName = names
}
indices = match(commName, names)

sadType = c("meteSAD", "empirSAD")
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
        if ( server == 'unc') 
          system(paste('bsub -q week -M', memory, '-J', names[i],
                       '-o', log_file, 'Rscript spat_analysis.R',
                       S[i], N[i], 200, bisect_fine[i], bisect_coarse[i],
                       grain_fine[i], FALSE, k, m, bisect, NA, NA,
                       name, TRUE, iteration, sep=' '))
        else
          system(paste('nice Rscript spat_analysis.R', S[i], N[i], 200,
                       bisect_fine[i], bisect_coarse[i],
                       grain_fine[i], FALSE, k, m, bisect, NA, NA,
                       name, TRUE, iteration, '>', log_file, '2>&1 &', sep=' '),
                 wait=FALSE)
      }    
    }
  }
}
