## Purpose: to compute the empirical observed DDR 

setwd('~/maxent/spat/scripts')

clArgs = commandArgs(trailingOnly=TRUE)

server = clArgs[1]
npar = clArgs[2]
nperm = clArgs[3]
commName = clArgs[4]

metricsToCalc = "sorensen"

dataType = c("abu", "binary")

for (i in commName) {
  for (j in metricsToCalc) {
    for (k in dataType) {
      log_file = paste('./log_files/error_', j, '_', k, '_', i,
                       '.log', sep='')
      if (server == 'unc') {
        if (k == 'abu') {
          system(paste('export OMP_NUM_THREADS=', npar, sep=''))
          system(paste('bsub -q week -M 20 -n', npar, "-R 'span[hosts=1]' -J",
                       i, '-o' log_file, 'Rscript spat_empir_analysis.R',
                       i, j, k, nperm, npar, sep=' '))
        }
        else {
          system('export OMP_NUM_THREADS=1')
          system(paste('bsub -q week -M 20 -n 1 -J', i, '-o' log_file,
                       'Rscript spat_empir_analysis.R', i, j, k,
                       'NULL', 1, sep=' '))
        }
      }
      else {
        if (k == 'abu')
          system(paste('Rscript spat_empir_analysis.R', i, j, k,
                       nperm, npar, '>' log_file, '2>&1 &', sep=' '),
                 wait=FALSE)
        else 
          system(paste('Rscript spat_empir_analysis.R', i, j, k,
                       'NULL', 1, '>' log_file, '2>&1 &', sep=' '),
                 wait=FALSE)
      }
    }
  }
}

