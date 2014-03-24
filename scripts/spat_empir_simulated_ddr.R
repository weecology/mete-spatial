## Purpose: to compute the empirical simulated DDR 

print('Submitting jobs to compute the METE simulated DDRs, ...')

clArgs = commandArgs(trailingOnly=TRUE)

server = clArgs[1]
commName = clArgs[2]
ncomm = clArgs[3]
dataType = clArgs[4]
bisect = clArgs[5]
memory = as.numeric(clArgs[6])
iteration = as.numeric(clArgs[7])

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
S = read.table('./data/S_vals.txt')
N = read.table('./data/N_vals.txt')
bisect_fine = read.table('./data/bisect_fine.txt')
bisect_coarse = read.table('./data/bisect_coarse.txt')
grain_fine = read.table('./data/grain_fine.txt')
names = as.character(read.table('./data/shrtnames.txt', colClasses='character'))

commName = unlist(strsplit(commName, ' '))
if (commName[1] == 'all') {
  commName = names
}
indices = match(commName, names)

sadType = c("meteSAD", "empirSAD")
metrics = "sorensen"

nice = FALSE # make TRUE to use nice when submitting jobs
wait = TRUE  # make FALSE to submit all the jobs at once

for (i in indices) {
  for (j in sadType) {
    for (k in dataType) {
      for (m in metrics) {
        log_file = paste('./scripts/log_files/error_sim_analysis_', names[i],
                         '_', j, '.log', sep='')
        if (j == 'meteSAD')
          name = names[i]
        if (j == 'empirSAD')
          name = paste(names[i], '_empirSAD', sep='')
        if (server == 'LSF') {
          system(paste('bsub -q week -M', memory, '-J', names[i],
                       '-o', log_file, 'Rscript ./scripts/spat_analysis.R',
                       S[i], N[i], ncomm, bisect_fine[i], bisect_coarse[i],
                       grain_fine[i], FALSE, k, m, bisect, NA, NA,
                       name, TRUE, iteration, sep=' '))
        }
        else {
          if (nice)
            cmd = 'nice Rscript'
          else
            cmd = 'Rscript'
          system(paste(cmd, './scripts/spat_analysis.R', S[i], N[i], ncomm,
                       bisect_fine[i], bisect_coarse[i],
                       grain_fine[i], FALSE, k, m, bisect, NA, NA,
                       name, TRUE, iteration, '>', log_file, '2>&1', sep=' '),
                 wait=wait)
        }  
      }    
    }
  }
}

print('Submitting jobs to compute the METE simulated DDRs, complete!')
