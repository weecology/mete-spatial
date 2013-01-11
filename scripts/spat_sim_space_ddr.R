setwd('~/maxent/spat/scripts')

clArgs = commandArgs(trailingOnly=TRUE)
server = clArgs[1] ## unc or usu

Svals = round(10^seq(log10(10), log10(100), length.out=20))
Nvals = round(10^seq(log10(120), log10(5e5), length.out=20))

ncomm = 200
bisec_fine = 12
bisec_coarse = 4
grain_fine = 1
transect = FALSE
dataType = 'abu'
metricsToCalc = 'sorensen'
direction = NA
tolerance = NA
name = NA
big = TRUE

for (S in Svals) {
  for (N in Nvals) {
    for (type in dataType) {
      log_file = paste('./log_files/error_', S, '_', N, '.log', sep='')
      if (server == unc)
        system(paste('bsub -q week -M 8 -J space -o', log_file,
                     'Rscript spat_analysis.R', S, N,
                     ncomm, bisec_fine, bisec_coarse, grain_fine, transect,
                     type, metricsToCalc, direction, tolerance, name, big,
                     sep=' '), wait = FALSE)
      if (server == usu)0
        system(paste('Rscript spat_analysis.R', S, N, ncomm, bisec_fine,
                     bisec_coarse, grain_fine, transect, type, 
                     metricsToCalc, direction, tolerance, name, big,
                     '>', log_file, '2>&1 &', sep=' '), wait = FALSE)
  }
}

