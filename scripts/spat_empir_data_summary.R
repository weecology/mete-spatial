## Author: Dan McGlinn
## Purpose: to generate a table that summarizes characteristics of the empirical
## datasets that are employed in the distance decay study

print('Exporting empirical data summaries, ...')

shrtnames = as.character(as.matrix(read.table('./data/shrtnames.txt')))

grain_fine = grain_coarse = bisect_fine = bisect_coarse = S = N = shape = NULL

comms = vector("list", length=length(shrtnames))
names(comms) = shrtnames
for (i in seq_along(comms)) {
  comms[[i]] = read.csv(paste('./data/', shrtnames[i], '_comms.csv', sep=''))
  grains = unique(comms[[i]][ , 1])
  grain_fine[i] = grains[1]
  grain_coarse[i] = tail(grains, 1)
  bisect_fine[i] = log2(nrow(comms[[i]][ comms[[i]][ , 1] == grain_fine[i], ]))
  bisect_coarse[i] = log2(nrow(comms[[i]][ comms[[i]][ , 1] == grain_coarse[i], ]))
  sad = read.csv(paste('./data/', shrtnames[i], '_sad.csv', sep=''), header=F)
  S[i] = ncol(sad)
  N[i] = sum(sad)
  shape[i] = ifelse(bisect_coarse[i] %% 2 == 0 , 'square', 'rectangle')
}

dat = data.frame(datnames = c('BCI', rep('Cocoli', 2), 'Crosstimbers',
                              rep('Sherman', 3),'Serpentine', 'Oosting',
                              'UCSC', 'Luquillo', 'Graveyard','Landsend',
                              'Rocky', 'Bormann', 'Wood Bridge', 'Bald Mnt', 
                              'Bryan', 'Big Oak'),
                 shrtname = c('bci', paste('cocoli', 1:2, sep=''), 'cross',
                              paste('sherman', 1:3, sep=''), 'serp', 'oosting',
                              'ucsc', 'luquillo', 'graveyard', 'landsend',
                              'rocky', 'bormann', 'woodbridge', 'baldmnt',
                              'bryan', 'bigoak'),
                 habitat = c(rep('tropical', 3), 'oak savanna',
                             rep('tropical', 3), 'grassland',  'oak-hickory',
                             'mixed evergreen', 'tropical', rep('pine',2), 
                             rep('oak-hickory', 6)),
                 area = c(50 * 1e4, rep(200 * 100, 2), 4 * 1e4,
                          rep(200 * 100, 2), 140^2, 64, 256^2, 150 * 300,
                          250 * 500, 100^2, 130 * 65, 120^2, 140^2, 71^2,
                          50 * 100, 185 * 92.5, 200 * 100) * 1e-4)

dat = dat[match(shrtnames, dat$shrtname), ]

grain_fine = dat$area / 2^bisect_fine * 1e4
grain_coarse = dat$area / 2^bisect_coarse * 1e4

## output overall summary
datSummary = data.frame(dat, shape,
                        bisect_fine, bisect_coarse, 
                        grain_fine, grain_coarse,
                        S, N)
names(datSummary) = c('site', 'shrtnames', 'habitat', 'area_ha', 'shape',
                      'bisect_fine', 'bisect_coarse', 'grain_fine_m2', 
                      'grain_coarse_m2', 'S', 'N')
write.csv(datSummary, file='empir_data_summary.csv', row.names=FALSE)

## output S, N, and bisection info 
write.table(matrix(S, nrow=1), file=file.path('./data', 'S_vals.txt'), sep=' ', 
            row.names=FALSE, col.names=FALSE)
write.table(matrix(N, nrow=1), file=file.path('./data', 'N_vals.txt'), sep=' ', 
            row.names=FALSE, col.names=FALSE)
write.table(matrix(bisect_fine, nrow=1), file=file.path('./data', 'bisect_fine.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE)
write.table(matrix(bisect_coarse, nrow=1), file=file.path('./data', 'bisect_coarse.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE)
write.table(matrix(grain_fine, nrow=1), file=file.path('./data', 'grain_fine.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE)
write.table(matrix(dat$habitat, nrow=1), file=file.path('./data', 'habitat.txt'),
            sep=' ', row.names=FALSE, col.names=FALSE)

print('Exporting empirical data summaries, complete!')

