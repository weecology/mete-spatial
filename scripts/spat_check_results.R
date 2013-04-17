

setwd('./sorensen')

files = dir()

sites = as.character(read.table('../data/shrtnames.txt', colClasses='character'))

load_list = function(filenames, R_obj=FALSE, obj_name=NULL) {
  dat = vector('list', length(filenames))
  for (i in seq_along(filenames)) {
    if (R_obj) {
      load(filenames[i])
      dat[[i]] = eval(parse(text=obj_name))
    }
    else {
      dat[[i]] = read.csv(filenames[i])
    }
  }
  return(dat)
}


## bisected binary empirical results

tmp = files[grep('_bisect_binary.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp, invert=TRUE)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites]

## bisected abundance empirical results

tmp = files[grep('_bisect_abu.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp, invert=TRUE)]
tmp = tmp[grep('_uni_', tmp, invert=TRUE)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites] 
## sherman1, sherman2, luquillo, ferp on killdevil
## landsend on jayne 
## bryan on wash
## bigoak on zoe

dat = load_list(tmp, TRUE, 'sorensen')
for(i in seq_along(dat)) { 
  print(paste(tmp[i], ncol(dat[[1]][[1]]$sorensenNull$vario)))
}


## bisected abundance univariate empirical results
tmp = files[grep('_bisect_abu.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp, invert=TRUE)]
tmp = tmp[grep('_uni_', tmp)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites] 

## mete analytical logseries results
tmp = files[grep('.csv', files)]
tmp = tmp[grep('_empirSAD_', tmp, invert=T)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][1]))
sites[!sites %in% tmp_sites]

## mete analytical empirSAD results
tmp = files[grep('.csv', files)]
tmp = tmp[grep('_empirSAD_', tmp)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][1]))
sites[!sites %in% tmp_sites]

## bisected binary logseries simulated results
tmp = files[grep('_bisect_binary.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp)]
tmp = tmp[grep('_empirSAD_', tmp, invert=T)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites] 

dat = load_list(tmp, TRUE, 'metrics')
for(i in seq_along(dat)) { 
  print(paste(tmp[i], sum(sapply(dat[[1]], function(x) !is.null(x)))))
}

## bisected abundance logseries simulated results
tmp = files[grep('_bisect_abu.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp)]
tmp = tmp[grep('_empirSAD_', tmp, invert=T)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites] 

## bisected binary empirSAD simulated results
tmp = files[grep('_bisect_binary.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp)]
tmp = tmp[grep('_empirSAD_', tmp)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites] 

## bisected abundance empirSAD simulated results
tmp = files[grep('_bisect_abu.Rdata', files)]
tmp = tmp[grep('_C200_B13_', tmp)]
tmp = tmp[grep('_empirSAD_', tmp)]
tmp_sites = as.character(sapply(tmp, function(x) strsplit(x, '_')[[1]][2]))
sites[!sites %in% tmp_sites]
