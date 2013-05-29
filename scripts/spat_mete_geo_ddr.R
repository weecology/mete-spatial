## Purpose: to compute the analytical mete DDR for geographic distance
## rather than seperation distance

setwd('~/maxent/spat')

source('./scripts/spat_functions.R')

shrtnames = read.table('./data/shrtnames.txt', colClasses='character')
grains = as.numeric(read.table('./data/grain_fine.txt'))

file_names = dir('./sorensen')
empirSAD_files = sort(file_names[grep('empirSAD_mete', file_names)])
sor_empirSAD = vector('list', length(empirSAD_files))
names(sor_empirSAD) = sub('_empirSAD_mete_sor.csv', '', empirSAD_files)
for (i in seq_along(sor_empirSAD))
  sor_empirSAD[[i]] = read.csv(file.path('./sorensen', empirSAD_files[i]))

for (i in seq_along(shrtnames)) {
  name = shrtnames[i]
  index = grep(name, names(sor_empirSAD))
  if (length(index) > 1) {
    sor_tmp = rbind(sor_empirSAD[[index[1]]], sor_empirSAD[[index[2]]])
    sor_tmp = sor_tmp[order(sor_tmp$i, sor_tmp$j, decreasing=T), ]
    write.csv(sor_tmp, 
              paste('./ddr/', name, '_empirSAD_mete_sor.csv', sep=''),
              row.names=F)
  }            
}

file_names = dir('./ddr')
empirSAD_files = sort(file_names[grep('empirSAD_mete', file_names)])
sor_empirSAD = vector('list', length(empirSAD_files))
names(sor_empirSAD) = sub('_empirSAD_mete_sor.csv', '', empirSAD_files)
for (i in seq_along(sor_empirSAD))
  sor_empirSAD[[i]] = read.csv(file.path('./ddr', empirSAD_files[i]))


## average and drop unneeded cocoli and sherman plots
sherman = c('sherman1', 'sherman2')
cocoli = c('cocoli1', 'cocoli2')
sor_empirSAD$sherman1$sor = apply(sapply(sor_empirSAD[sherman], 
                                         function(x) x$sor), 1, mean)
sor_empirSAD$cocoli1$sor = apply(sapply(sor_empirSAD[cocoli], 
                                         function(x) x$sor), 1, mean)
indices_to_drop = match(c('sherman2', 'sherman3', 'cocoli2'), names(sor_empirSAD))
sor_empirSAD = sor_empirSAD[-indices_to_drop]

m = 1
n = 4
bisections = unique(sor_empirSAD[[m]]$i)
i_bisect = bisections[n]
sep_dist = as.matrix(get_sep_dist(i_bisect))
sor_vals = sor_empirSAD[[m]]$sor[sor_empirSAD[[m]]$i == i_bisect]
sor_vals = sor_vals[order(sor_empirSAD[[m]]$j[sor_empirSAD[[m]]$i == i_bisect])]
sor_vals = c(1, sor_vals)
sor_dist = sor_vals[as.numeric(sep_dist + 1)]
sor_mat = matrix(sor_dist, 32, 32, byrow=FALSE)
sor_dist = as.dist(sor_mat)
coords = get_bisect_coords(i_bisect)
grain = grains[match(names(sor_empirSAD)[m], shrtnames)]
grain_mult = 2 ^ i_bisect
coords[ , 1:2] = coords[ , 1:2] * sqrt(grain * grain_mult)
geo_dist = dist(coords[ , 1:2])
geo_dist = ifelse(geo_dist > max(geo_dist) / 2, NA, geo_dist)

plot(sor ~ dist, data=sor_empirSAD[[m]], subset = i == 5, log='xy',
     type='l', ylim=c(.1,1), xlim=c(5, 80))
points(geo_dist, sor_dist)
lines(lowess(geo_dist[!is.na(geo_dist)],
             sor_dist[!is.na(geo_dist)]), col='red')
