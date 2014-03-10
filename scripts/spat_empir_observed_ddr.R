## Purpose: to compute the empirical observed DDR 

print('Submitting jobs to compute the empirical DDRs, ...')

clArgs = commandArgs(trailingOnly=TRUE)

server = clArgs[1] ##'unc' or something else
npar = clArgs[2]
nperm = clArgs[3]
commName = clArgs[4]
method = clArgs[5] ## 'multi' or 'uni' for multivarite and univariate respectively
dataType = clArgs[6] ## 'abu', 'binary', or 'both'

metricsToCalc = "sorensen"

if (dataType == 'both') {
  dataType = c("abu", "binary")
}

commName = unlist(strsplit(commName, ' '))
if (commName[1] == 'all') {
  commName = as.character(read.table('./data/shrtnames.txt',
                                     colClasses='character')[1, ])
}


for (i in commName) {
  for (j in metricsToCalc) {
    for (k in dataType) {
      log_file = paste('./scripts/log_files/error_', j, '_', k, '_', i,
                       '.log', sep='')
      if (server == 'unc') {
        if (k == 'abu') {
          system(paste('export OMP_NUM_THREADS=', npar, sep=''))
          system(paste('bsub -q week -M 40 -n', npar, "-R 'span[hosts=1]' -J",
                       i, '-o', log_file, 'Rscript ./scripts/spat_empir_analysis.R',
                       i, j, k, method, npar, nperm, sep=' '))
        }
        else {
          system('export OMP_NUM_THREADS=1')
          system(paste('bsub -q week -M 7 -n 1 -J', i, '-o', log_file,
                       'Rscript ./scripts/spat_empir_analysis.R', i, j, k,
                       method, 1, NA, sep=' '))
        }
      }
      else {
        if (k == 'abu')
          system(paste('nice Rscript ./scripts/spat_empir_analysis.R', i, j, k,
                       method, npar, nperm, '>', log_file, '2>&1',
                       sep=' '), wait=FALSE)
        else 
          system(paste('nice Rscript ./scripts/spat_empir_analysis.R', i, j, k,
                       method, 1, NA, '>', log_file, '2>&1', sep=' '),
                       wait=FALSE)
      }
    }
  }
}

print('Submitting jobs to compute the empirical DDRs, complete!')
