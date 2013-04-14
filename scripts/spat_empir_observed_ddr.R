## Purpose: to compute the empirical observed DDR 

setwd('~/maxent/spat/scripts')

clArgs = commandArgs(trailingOnly=TRUE)

server = clArgs[1]
npar = clArgs[2]
nperm = clArgs[3]
commName = clArgs[4]
method = clArgs[5]
dataType = clArgs[6]

metricsToCalc = "sorensen"

if (dataType == 'both') {
  dataType = c("abu", "binary")
}

if (commName == 'all') {
  commName = as.character(read.table('../data/shrtnames.txt',
                                     colClasses='character')[1,])
}

for (i in commName) {
  for (j in metricsToCalc) {
    for (k in dataType) {
      log_file = paste('./log_files/error_', j, '_', k, '_', i,
                       '.log', sep='')
      if (server == 'unc') {
        if (k == 'abu') {
          system(paste('export OMP_NUM_THREADS=', npar, sep=''))
          system(paste('bsub -q week -M 40 -n', npar, "-R 'span[hosts=1]' -J",
                       i, '-o', log_file, 'Rscript spat_empir_analysis.R',
                       i, j, k, method, npar, nperm, sep=' '))
        }
        else {
          system('export OMP_NUM_THREADS=1')
          system(paste('bsub -q week -M 7 -n 1 -J', i, '-o', log_file,
                       'Rscript spat_empir_analysis.R', i, j, k,
                       method, 1, NA, sep=' '))
        }
      }
      else {
        if (k == 'abu')
          system(paste('nice Rscript spat_empir_analysis.R', i, j, k,
                       method, npar, nperm, '>', log_file, '2>&1',
                       sep=' '), wait=FALSE)
        else 
          system(paste('nice Rscript spat_empir_analysis.R', i, j, k,
                       method, 1, NA, '>', log_file, '2>&1', sep=' '),
                       wait=FALSE)
      }
    }
  }
}

