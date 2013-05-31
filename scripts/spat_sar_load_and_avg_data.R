## Loads and averages the species area anaylsis data products

setwd('~/maxent/spat')

source('./scripts/spat_functions.R')

print('Aggregating SAR results, ...')

fileNames = dir('./sar')

## empirical SARs
empirFiles = grep('empir_sar.csv', fileNames, value=TRUE)
empir = vector('list', length(empirFiles))
names(empir) = sub('_empir_sar.csv', '', empirFiles)
for (i in seq_along(empirFiles))
  empir[[i]] = read.csv(paste('./sar/', empirFiles[i], sep=''))

## mete SARs analytical
meteFiles = grep('mete_sar.txt', fileNames, value=TRUE)
meteEmpirSADFiles = grep('empirSAD', meteFiles, value=TRUE)
meteLogSerFiles = grep('empirSAD', meteFiles, value=TRUE, invert=TRUE)
meteLogSer = vector('list', length(meteLogSerFiles))
meteEmpirSAD = vector('list', length(meteEmpirSADFiles))
names(meteLogSer) = sub('_mete_sar.txt', '', meteLogSerFiles)
names(meteEmpirSAD) = sub('_empirSAD_mete_sar.txt', '', meteEmpirSADFiles)
for (i in seq_along(meteLogSerFiles)) {
  meteLogSer[[i]] = read.csv(paste('./sar/', meteLogSerFiles[i], sep=''))
  meteLogSer[[i]]$area = meteLogSer[[i]]$area * empir[[i]]$area[1]
}
for (i in seq_along(meteEmpirSADFiles)) {
  meteEmpirSAD[[i]] = read.csv(paste('./sar/', meteEmpirSADFiles[i], sep=''))
  meteEmpirSAD[[i]]$area = meteEmpirSAD[[i]]$area * empir[[i]]$area[1]
}

## load expected sars under random placement, object is called 'srExp'
load('./sar/expected_empir_sars.Rdata')

sites = names(empir)
if (all(as.vector(sapply(c('sherman', 'cocoli'), grepl, sites)))) {
  ## average cocoli and sherman plots
  sherman = c('sherman1', 'sherman2')
  cocoli = c('cocoli1', 'cocoli2')
  empir = avg_site_results(empir, sherman)
  empir = avg_site_results(empir, cocoli)
  empir = empir[-match('sherman3', names(empir))]
  
  meteLogSer = avg_site_results(meteLogSer, sherman)
  meteLogSer = avg_site_results(meteLogSer, cocoli)
  meteLogSer = meteLogSer[-match('sherman3', names(meteLogSer))]
  
  meteEmpirSAD = avg_site_results(meteEmpirSAD, sherman)
  meteEmpirSAD = avg_site_results(meteEmpirSAD, cocoli)
  meteEmpirSAD = meteEmpirSAD[-match('sherman3', names(meteEmpirSAD))]
  
  meteAvgLogSer = avg_site_results(meteAvgLogSer, sherman)
  meteAvgLogSer = avg_site_results(meteAvgLogSer, cocoli)
  meteAvgLogSer = meteAvgLogSer[-match('sherman3', names(meteAvgLogSer))]
  
  meteAvgEmpirSAD = avg_site_results(meteAvgEmpirSAD, sherman)
  meteAvgEmpirSAD = avg_site_results(meteAvgEmpirSAD, cocoli)
  meteAvgEmpirSAD = meteAvgEmpirSAD[-match('sherman3', names(meteAvgEmpirSAD))]
  
  srExp = avg_site_results(srExp, sherman)
  srExp = avg_site_results(srExp, cocoli)
  srExp = srExp[-match('sherman3', names(srExp))]
}

## compute residuals
meteLogSer_iter_res = get_sar_resids(empir, meteLogSer, 'richness', 'sr_iter')
meteLogSer_noniter_res = get_sar_resids(empir, meteLogSer, 'richness', 'sr_noniter')

meteEmpirSAD_iter_res = get_sar_resids(empir, meteEmpirSAD, 'richness', 'sr_iter')
meteEmpirSAD_noniter_res = get_sar_resids(empir, meteEmpirSAD, 'richness', 'sr_noniter')

rpLogSer_res = get_sar_resids(empir, srExp, 'richness', 'S_logser_binom')
rpEmpirSAD_res = get_sar_resids(empir, srExp, 'richness', 'S_binom')


sar_res = data.frame(site = meteLogSer_iter_res$site, area = meteLogSer_iter_res$area,
                     logser_iter = meteLogSer_iter_res$res,
                     logser_noniter = meteLogSer_noniter_res$res,
                     logser_rp = rpLogSer_res$res,
                     empirsad_iter = meteEmpirSAD_iter_res$res,
                     empirsad_noniter = meteEmpirSAD_noniter_res$res,
                     empirsad_rp = rpEmpirSAD_res$res)

## Add richness and individual density to residual data.frame
sar_data = read.csv('./sar/empir_sars.csv')
sar_res = merge(sar_res, 
                sar_data[ , c('site', 'area', 'richness', 'sr_std', 'indiv')],
                all.x=TRUE)
## sort sar_res so that its in a reasonable order
sar_res = sar_res[order(sar_res$site, sar_res$area), ]

## bring in habitat type
shrtnm = as.character(read.table('./data/shrtnames.txt', colClasses='character'))
habitat = as.character(read.table('./data/habitat.txt', colClasses='character'))
sar_res$hab = habitat[match(sar_res$site, shrtnm)]

sar_data = sar_res
sar_data[ , 3:8] = sar_data$richness - sar_data[ , 3:8] 

## export files
write.csv(sar_res, file='./sar/sar_residuals.csv', row.names=F)
write.csv(sar_data, file='./sar/sar_raw_data.csv', row.names=F)

print('Aggregating SAR results, complete!')
