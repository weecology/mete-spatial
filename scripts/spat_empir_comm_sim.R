## Purpose: to simulate communities that are parameterized by the empirically 
## observed communities.

setwd('~/maxent/spat/scripts')

clArgs = commandArgs(trailingOnly=TRUE)
server = clArgs[1]
indices=clArgs[2]

names = read.table('../data/shrtnames.txt', colClasses='character')
bisect = read.table('../data/bisect_fine.txt')
S = read.table('../data/S_vals.txt') 
N = read.table('../data/N_vals.txt') 

sadType = c("meteSAD", "empirSAD")

for (i in indices) {
  for (j in sadType) {
    log_file = paste('./log_files/error_sim_analysis_', names[i],
                     '_', j, '.log', sep='')
    if (j == 'meteSAD')
      abu_file = 'None'
    else
      abu_file = paste('../data/', names[i], '_sad.csv', sep='')
    if (server == 'unc')
      system(paste('bsub -q week -M 8 -J', names[i], '-o',
                   log_file, 'python spat_community_generation.py',
                   S[i], N[i], 200, bisect[i], 'False', 'None',
                   names[i], sep=' '))
    else
      system(paste('nice python spat_community_generation.py',
                   S[i], N[i], 200, bisect[i], 'False', 'None',
                   names[i], '>', log_file, '2>&1', sep=' '),
             wait=FALSE)
  }
}
