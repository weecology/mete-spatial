
setwd('~/maxent/spat')
source('./scripts/spat_functions.R')

print('Computing binomial expected SARs, ...')

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
S0 = sapply(dat, length)
N0 = sapply(dat, sum)
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
  S_binom = sapply(grains[[i]], function(g) 
                   exp_S_binom(g, maxA, dat[[i]]))
  sd_S_binom = sapply(grains[[i]], function(g) 
                      sqrt(exp_Svar_binom(g, maxA, dat[[i]]))) 
  S_lo = S_binom - sd_S_binom
  S_hi = S_binom + sd_S_binom
  if (names(srExp)[i] != 'cross')
    S_logser_binom = sar_logser_binom(grains[[i]], maxA, S0[[i]], N0[[i]])
  else
    S_logser_binom = NA
  srExp[[i]]  = data.frame(grains = grains[[i]], S_binom, S_lo, S_hi, S_logser_binom)
}

save(srExp, file='./sar/expected_empir_sars.Rdata')

print('Computing binomial expected SARs, complete!')
