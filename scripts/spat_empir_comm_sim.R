## Purpose: to simulate communities that are parameterized by the empirically 
## observed communities.

print('Submitting jobs to simulate empirical communities under METE, ...')

clArgs = commandArgs(trailingOnly=TRUE)
server = clArgs[1]
ncomm = clArgs[2]
sites = clArgs[-(1:2)]

names = read.table('./data/shrtnames.txt', colClasses='character')
bisect = read.table('./data/bisect_fine.txt')
S = read.table('./data/S_vals.txt') 
N = read.table('./data/N_vals.txt') 

sadType = c("meteSAD", "empirSAD")

indices = match(sites, names)

nice = FALSE # make TRUE to use nice when submitting jobs
wait = TRUE  # make FALSE to submit all the jobs at once

for (i in indices) {
  for (j in sadType) {
    log_file = paste('./scripts/log_files/error_mete_comm_gen_', names[i],
                     '_', j, '.log', sep='')
    if (j == 'meteSAD')
      abu_file = 'None'
    else
      abu_file = paste('./data/', names[i], '_sad.csv', sep='')
    if (server == 'LSF') {
      system(paste('bsub -q week -M 8 -J', names[i], '-o',
                   log_file, 'python ./scripts/spat_community_generation.py',
                   S[i], N[i], ncomm, bisect[i], 'False', abu_file,
                   names[i], sep=' '))
    }
    else {
      if (nice)
        cmd = 'nice python'
      else
        cmd = 'python'
      system(paste(cmd, './scripts/spat_community_generation.py',
                   S[i], N[i], ncomm, bisect[i], 'False', abu_file,
                   names[i], '>', log_file, '2>&1', sep=' '),
             wait=wait)
    }  
  }
}

print('Submitting jobs to simulate empirical communities under METE, complete!')

