
setwd('~/maxent/spat')

source('./scripts/spat_sim_vario_func.R')

shrtnames = c('bci', 'cocoli1', 'cocoli2', 'cross', 'sherman1', 'sherman2', 
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak')

dat = vector("list", length=length(shrtnames))
names(dat) = shrtnames
for (i in seq_along(dat))
  dat[[i]] = read.csv(paste('./data/', shrtnames[i], '_sad.csv', sep=''),
                      header=FALSE)

shrtnm = read.table('./data/shrtnames.txt', colClasses='character')
bisect = read.table('./data/bisect_fine.txt')
grain_fine = read.table('./data/grain_fine.txt')
grains = vector("list", length=length(shrtnames))
names(grains) = shrtnames
for (i in seq_along(grains)) {
  index = match(names(dat[i]), shrtnm)
  B = as.numeric(bisect[index])
  Amin = as.numeric(grain_fine[index])
  grains[[i]] = 2^(0:B) * Amin
}

srExp = vector("list", length=length(grains))
names(srExp) = shrtnames
for (i in seq_along(srExp)) {
  maxA = max(grains[[i]])
  srCol = sapply(grains[[i]], function(g) 
                 spAvgExpColeman(g / maxA, dat[[i]]))
  sdCol = sapply(grains[[i]], function(g) 
                 sqrt(spVarExpColeman(g / maxA, dat[[i]]))) 
  sr.lo = srCol - sdCol
  sr.hi = srCol + sdCol
  srExp[[i]]  = data.frame(grains = grains[[i]], srCol, sr.lo, sr.hi)
}

save(srExp, file='./sar/expected_empir_sars.Rdata')

plot(srCol ~ grains, data=srExp[[1]], type='o', log='xy')
lines(sr.lo ~ grains, data=srExp[[1]], col='red')
lines(sr.hi ~ grains, data=srExp[[1]], col='red')



  